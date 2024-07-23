library(Seurat)
library(DoubletFinder)
library(dplyr)
library(magrittr)
library(ggplot2)
############################## Read raw data ###################################
setwd("counts/")
folders=list.files('.')
Raw = lapply(folders,function(folder){CreateSeuratObject(counts = Read10X_h5(folder), project = substr(folder,1,3),min.cells = 0,min.features = 0)})
saveRDS(Raw,file = 'Raw.rds')
############################### Filter data ####################################
NormList <- Raw
npcs = 30
for (i in 1:length(NormList)) {
  Object <- NormList[[i]]

  ### filter low/high gene/counts cell
  qt_Feature = quantile(Object$nFeature_RNA,c(0.025,0.975))
  qt_Counts = quantile(Object$nCount_RNA,c(0.025,0.975))
  Object <- subset(Object,subset = (nFeature_RNA>qt_Feature[1]&nFeature_RNA < qt_Feature[2]) &
                                   (nCount_RNA>qt_Counts[1]&nCount_RNA<qt_Counts[2]))

  ### filter cells with high expression of mitochondria 
  Object$percent.mt <- PercentageFeatureSet(Object, pattern = "^mt-")
  Object <- subset(Object,percent.mt<30)
  
  ### remove ribsome/mitochondria genes
  Object <- Object[!grepl("^Rp[ls]",rownames(Object)),]
  Object <- Object[!grepl("^mt",rownames(Object)),]
  
  Object$Condition <- substr(Object$orig.ident,1,2)
  Object$cellnames <- colnames(Object)
  
  Object <- NormalizeData(Object) %>%
            ScaleData() %>%
            FindVariableFeatures() %>%
            RunPCA(npcs = npcs)
  
  ### find doublets
  pK_bcmvn <- paramSweep_v3(Object, PCs = 1:npcs, sct = FALSE) %>%
              summarizeSweep(GT = FALSE) %>%
              find.pK() %$%
              pK[which.max(BCmetric)] %>%
              as.character() %>%
              as.numeric()
  nExp_poi <- round(0.075*nrow(Object@meta.data))
  Object <- doubletFinder_v3(Object,
                                    PCs = 1:npcs,
                                    pN = 0.25,
                                    pK = pK_bcmvn,
                                    nExp = nExp_poi,
                                    reuse.pANN = FALSE,
                                    sct = FALSE)
  names(Object@meta.data)[ncol(Object@meta.data)] <- "Doublets"
  Object <- subset(Object,Doublets == 'Singlet')
  NormList[[i]] <- Object
}
sum(sapply(NormList,dim)[2,])
names(NormList) <- c("KO1","KO2","KO3","WT1","WT2","WT3")
############################# Integrate Data ###################################
### Integrate
features <- SelectIntegrationFeatures(object.list = NormList)
Anchors <- FindIntegrationAnchors(object.list = NormList,anchor.features = features,dims = 1:npcs)
Integ<- IntegrateData(anchorset = Anchors, dims = 1:npcs) 
DefaultAssay(Integd) <- "integrated"

### Find variable genes from expression level bins
source('./FindVarBins.R')
Integ$project <- 'Integrated'
Integ <- FindVarBins(Integ,assay.value = 'RNA')
Integ <- ScaleData(Integ, verbose = FALSE) 
npcs=30
Integ <- RunPCA(Integ,npcs = npcs, verbose = FALSE) %>%
          RunUMAP(dims = 1:npcs,reduction = 'pca') %>%
          FindNeighbors(dims = 1:npcs, reduction="pca")
ElbowPlot(Integ,ndims = npcs)
DimPlot(Integ,reduction = 'umap',group.by = "Condition")
########################## Annotate cell type ##################################
sce <- Integ
DefaultAssay(sce)<-'integrated'
sce$Condition <- factor(sce$Condition,levels = c('WT','KO'))
sce$orig.ident <- factor(sce$orig.ident,levels = c('WT1','WT2','WT3','KO1','KO2','KO3'))
sce <- FindClusters(sce,algorithm = 2,resolution = 0.7)
DimPlot(sce,reduction = 'umap',label = T,repel = T,label.box = T,group.by = 'seurat_clusters')+NoLegend()

DefaultAssay(sce)<-'RNA'
signatures = read.csv("./Markers.csv",header = T,na.strings = "") %>%
  as.list() %>%
  lapply(function(x){x <- x[!is.na(x) & x %in% rownames(sce@assays$RNA)]})
signatures = signatures[lengths(signatures)!=0]
signatures
names(signatures) <- gsub('\\.','-',names(signatures))

DotPlot(sce,features = signatures,scale = T,assay = "RNA")+
  scale_color_gradient2(low = "white",high = "blue",midpoint = 0) +
  RotatedAxis()

celltype <- as.character(sce$seurat_clusters)
celltype[celltype %in% c(22)] <- 'UE'
celltype[celltype %in% c(21)] <- 'Podo'
celltype[celltype %in% c(9)] <- 'ALH'
celltype[celltype %in% c(3,24)] <- 'DLH'
celltype[celltype %in% c(13)] <- 'DCT'
celltype[celltype %in% c(18)] <- 'CD-PC'
celltype[celltype %in% c(17)] <- 'CD-IC'
celltype[celltype %in% c(8,10,16)] <- 'Endo'
celltype[celltype %in% c(23,15)] <- 'Stromal'
celltype[celltype %in% c(25)] <- 'RBC' 
celltype[celltype %in% c(0:2,4:7,11,12,14,19,20)] <- 'PT' 

sce$celltype <- factor(celltype,levels = c('PT','DCT','DLH','ALH','CD-IC','CD-PC','Podo','Endo','Stromal','UE','RBC'))
Idents(sce) <- 'celltype'
DimPlot(sce,label = T)

DEGs <- FindAllMarkers(sce,only.pos = T,return.thresh = 0.05,logfc.threshold = 0.5,min.diff.pct = 0.2)

############################ Recluster PT ###################################
PT <- subset(sce,celltype %in% c('PT'))
DefaultAssay(PT) <- 'integrated'
PT <- FindVarBins(PT,ntop = 25) %>% 
      ScaleData() %>% 
      RunPCA(npcs=20)
  
PT <- RunUMAP(PT,dims=1:20) %>%
      FindNeighbors(dims = 1:20)
PT <- FindClusters(PT,algorithm = 2,resolution = 0.3)

DimPlot(PT,label = T)

DefaultAssay(PT) <- 'RNA'
PT@assays$RNA@var.features <- PT@assays$integrated@var.features

### transfer ref
ref <- readRDS('GSE151658_integrated.0h.1h_4h_16h_27h_36h_48h.rds')
celltype <- as.character(ref$integrated_snn_res.1)
celltype[celltype %in% c(0,3,12)] <- 'S1'
celltype[celltype %in% c(1,8,11)] <- 'S3'
celltype[celltype %in% c(2,6,7)] <- 'S2'
celltype[celltype %in% c(4)] <- 'DCT'
celltype[celltype %in% c(5)] <- 'LOH'
celltype[celltype %in% c(9)] <- 'Endo'
celltype[celltype %in% c(10)] <- 'Mφ/DC'
celltype[celltype %in% c(13)] <- 'Lymphocytes'
celltype[celltype %in% c(14)] <- 'CNT'
celltype[celltype %in% c(15)] <- 'CD-PC'
celltype[celltype %in% c(16)] <- 'CD-IC'
celltype[celltype %in% c(17)] <- 'S3T2'
celltype[celltype %in% c(18,22,24,25)] <- 'Other'
celltype[celltype %in% c(19)] <- 'Neutrophil'
celltype[celltype %in% c(20)] <- 'Peri/Stromal'
celltype[celltype %in% c(21)] <- 'Pro'
celltype[celltype %in% c(23)] <- 'Other/UE'
ref$celltype <- factor(celltype)
DimPlot(ref,group.by = 'celltype',label = T)

ref_use <- subset(ref,hrs=='0hr' & celltype %in% c('S1','S2','S3','S3T2'))
DimPlot(ref_use,group.by = 'celltype',label = T)

ref_use$project <- 'ref'
ref_use <- FindVarBins(ref_use,ntop = 25)
DefaultAssay(ref_use) <- 'RNA'
ref_use@assays$RNA@var.features <- ref_use@assays$integrated@var.features

anchors <- FindTransferAnchors(reference = ref_use, query = PT,query.assay = "RNA",reference.assay='RNA')
predictions <- TransferData(anchorset = anchors, refdata = ref_use$celltype)
celltype_trans <- predictions$predicted.id
names(celltype_trans) <- rownames(predictions)

PT$celltype_pre <- celltype_trans
DimPlot(PT,group.by = 'celltype_pre',label = T)

table(PT$celltype_pre,PT$seurat_clusters)

df1 = AverageExpression(ref_use,group.by = 'celltype',assays = 'RNA',slot = 'data',features = ref_use@assays$RNA@var.features)$RNA
colnames(df1) <- paste('ref',colnames(df1),sep='_')
df1
df2 = AverageExpression(PT,group.by = 'celltype_pre',assays = 'RNA',slot = 'data',features = PT@assays$RNA@var.features)$RNA
genes = intersect(rownames(df1),rownames(df2))
df2
df = cbind(df1[genes,],df2[genes,])

cor_mat = cor(df,method = 'spearman')[-(1:4),1:4]
cor_mat

### annotated subclsuter
celltype <- as.character(PT$seurat_clusters)
celltype[celltype %in% c(4,7)] <- 'PT-S1'
celltype[celltype %in% c(0,2,8)] <- 'PT-S2'
celltype[celltype %in% c(1,3,5)] <- 'PT-S3'
celltype[celltype %in% c(6)] <- 'PT-S3T2'

PT$celltype <- factor(celltype,levels = c('PT-S1','PT-S2','PT-S3','PT-S3T2'))
Idents(PT) <- 'celltype'
DimPlot(PT,label = T)
############################ Recluster Stromal #################################
library(reticulate)
library(sceasy)
library(Seurat)
use_python("/home/wangjing/miniconda3/envs/scvi-env/bin/python")
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
set.seed(123)

signatures_ST <- list('Peri'=c('Pdgfrb','Rgs5','Angpt2','Tagln','Cspg4'),
                      'Myo'=c('Acta2',"Col1a2","Plac8","Fn1"),
                      'Fibro'=c('Dcn',"Col1a1","Gsn","Fbln1","Lum","Cd34"))

Stromal <- subset(sce,celltype %in% c('Stromal'))
DefaultAssay(Stromal) <- 'RNA'
Stromal <- SCTransform(Stromal,vars.to.regress = c('nCount_RNA','nFeature_RNA','percent.mt'),)

Stromal <- CreateSeuratObject(counts = Stromal@assays$SCT@data,meta.data = Stromal@meta.data)
adata <- convertFormat(Stromal, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)

pcs = 10
scvi$settings$seed = as.integer(123)
scvi$model$SCVI$setup_anndata(adata,batch_key = "orig.ident")
model = scvi$model$SCVI(adata,n_latent=as.integer(pcs))
scvi$model$SCVI$view_anndata_setup(model)
model$train(early_stopping = T,early_stopping_monitor = 'train_loss')

latent = model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(Stromal)
Stromal[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(Stromal))
Stromal <- RunUMAP(Stromal, dims = 1:pcs, reduction = 'scvi')
Stromal <- FindNeighbors(Stromal, dims = 1:pcs, reduction = 'scvi')
Stromal <- FindClusters(Stromal,algorithm = 2,resolution = 0.2)
DimPlot(Stromal,label = T)

### annotated subclsuter
celltype <- as.character(Stromal$seurat_clusters)
celltype[celltype == 0] <- 'Unknow'
celltype[celltype == 1] <- 'Peri'
celltype[celltype == 2] <- 'Myo'
celltype[celltype == 3] <- 'Fibro'
Stromal$celltype <- factor(celltype,levels = c('Peri','Myo','Fibro','Unknow'))
Idents(Stromal) <- 'celltype'
DimPlot(Stromal,label = T)

######################### merge labels ################################
celltype_fine <- as.character(sce$celltype)
names(celltype_fine) <- colnames(sce)
celltype_fine[colnames(PT)] <- as.character(PT$celltype)
celltype_fine[colnames(Stromal)] <- as.character(Stromal$celltype)

sce$celltype_fine <- factor(celltype_fine,levels = c('PT-S1','PT-S2','PT-S3','PT-S3T2','DCT','DLH','ALH','CD-IC','CD-PC','Podo','Endo','Unknow','Peri','Fibro','Myo','UE','RBC'))
Idents(sce) <- 'celltype_fine'
DimPlot(sce,group.by = 'celltype_fine',label = T)
########################## senescence phenotype #######################
library(UCell)
DefaultAssay(sce) <- "RNA"
### Aging score
library(openxlsx)
AgeSigs <- read.xlsx("~/wangj2/GeneSets/msAgingScoreList.xlsx",sheet = 8)
AgeSigs <- AgeSigs$`Gene.(mouse)`
AgeSigs = AgeSigs[AgeSigs %in% rownames(sce)]
AgeSigs = list('Senescence' = AgeSigs)
sce <- AddModuleScore_UCell(sce,assay='RNA',slot = 'data',features=AgeSigs,ncores = 40,maxRank=2000)
################################ DEGs KO vs WT #######################################
Object <- sce
DefaultAssay(Object) <- "RNA"
Object<- ScaleData(Object,features = rownames(Object))

### Condition
Idents(Object) <- "Condition"
DEGs_Con <- FindMarkers(Object,only.pos = F,return.thresh = 0.05,assay = "RNA",ident.1 = "KO")

### Condition celltypes
Idents(Object) <- "celltype_fine"
clusters = levels(Object@active.ident)

DEGs_ConClus = lapply(clusters, function(x){
  print(x)
  if ( x %in% clusters) {
    DEG_ConClu <- FindMarkers(Object,ident.1 = "KO",group.by = "Condition",subset.ident = x,only.pos = F,
                              return.thresh = 0.05,assay = "RNA",logfc.threshold = 0,min.diff.pct = 0,min.pct = 0)
  }
})
names(DEGs_ConClus) <- clusters
sapply(DEGs_ConClus,dim)

################################ Enrichment ####################################
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

bg=bitr(rownames(sce@assays$RNA),'SYMBOL','ENTREZID','org.Mm.eg.db')

res = DEGs_Con
res = filter(DEGs_Con,avg_log2FC > 1 & p_val_adj < 0.05)
genelist <- bitr(rownames(res), fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)
ego <- enrichGO(gene = genelist, OrgDb = 'org.Mm.eg.db',universe = bg$ENTREZID,ont = "All")
go <- ego@result[,c("Description","qvalue","GeneRatio",'Count')]
go$Ratio <- as.numeric(unlist(lapply(strsplit(go$GeneRatio,"/"),"[",1)))/as.numeric(unlist(lapply(strsplit(go$GeneRatio,"/"),"[",2)))
df_up <- go %>% filter(qvalue < 0.05)


res = DEGs_Con
res = filter(res,avg_log2FC < -1 & p_val_adj < 0.05)
genelist <- bitr(rownames(res), fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)
ego <- enrichGO(gene = genelist, OrgDb = 'org.Mm.eg.db',universe = bg$ENTREZID,ont = "All")
go <- ego@result[,c("Description","qvalue","GeneRatio",'Count')]
go$Ratio <- as.numeric(unlist(lapply(strsplit(go$GeneRatio,"/"),"[",1)))/as.numeric(unlist(lapply(strsplit(go$GeneRatio,"/"),"[",2)))
df_down <- go %>% filter(qvalue < 0.05)
#################### ECM phenotype #########################
ECMSigs =c('Thbs2', 'Col4a4', 'Tnr', 'Col6a6', 'Vwf', 'Gp6', 'Col1a2', 'Itga3', 'Itga9', 'Col4a5', 'Gp5', 'Lamc1', 'Itga6', 'Itga1', 'Itgb5', 'Itgb8', 'Sv2c', 'Thbs3', 'Hmmr', 'Cd44', 'Lama3', 'Lamb1', 'Lamb3', 'Col6a4', 'Lama4', 'Itgb4', 'Lamc3', 'Itgb7', 'Col1a1', 'Fn1', 'Vtn', 'Sv2a', 'Col9a2', 'Itgb3', 'Gp1bb', 'Col9a3', 'Hspg2', 'Npnt', 'Gp1ba', 'Col9a1', 'Itga5', 'Itgav', 'Lama2', 'Col2a1', 'Spp1', 'Sdc4', 'Tnxb', 'Itga10', 'Col4a3', 'Sdc1', 'Lama5', 'Sv2b', 'Col6a3', 'Reln', 'Itga7', 'Itga2', 'Lamc2', 'Col6a2', 'Ibsp', 'Col4a6', 'Comp', 'Col6a1', 'Tnn', 'Itga11', 'Itgb6', 'Lamb2', 'Chad', 'Thbs4', 'Gp9', 'Itga8', 'Dag1', 'Itga2b', 'Agrn', 'Col4a2', 'Col4a1', 'Cd36', 'Lama1', 'Col6a5', 'Itgb1', 'Itga4', 'Thbs1', 'Tnc', 'Cd47')
ECMSigs = ECMSigs[ECMSigs %in% rownames(sce)]
ECMSigs = list('ECM' = ECMSigs)
sce <- AddModuleScore_UCell(sce,assay='RNA',slot = 'data',features=ECMSigs,ncores = 40,maxRank=2000)
############################## Cell type proportion #############################
library(ggpubr)
library(rstatix)

### all
freq = as.matrix(table(sce$orig.ident,sce$celltype))
freq = freq/rowSums(freq)
freq = data.frame(freq)
colnames(freq) <- c("Batch","Cluster","Percent")
freq$Condition <- factor(substr(freq$Batch,1,2))
head(freq)
pwc <- freq %>% group_by(Cluster) %>%
  t_test(Percent ~ Condition, paired = F)%>%
  add_significance("p")
pwc
pwc <- pwc %>% add_xy_position(x = "Cluster")
pwc$y.position = pwc$y.position-0.08 


### Stromal
freq = as.matrix(table(Stromal$orig.ident,as.character(Stromal$celltype)))
freq = freq/rowSums(freq)
freq = data.frame(freq)
colnames(freq) <- c("Batch","Cluster","Percent")
freq$Condition <- factor(substr(freq$Batch,1,2))
head(freq)
pwc <- freq %>% group_by(Cluster) %>%
  t_test(Percent ~ Condition, paired = F)%>%
  add_significance("p")
pwc
pwc <- pwc %>% add_xy_position(x = "Cluster")
############################## velocyto #############################
library(velocyto.R)
library(tidyverse)
library(SeuratWrappers)

velo <- read.loom.matrices("merged.loom")
velo <- as.Seurat(x = velo)

cellID_velo <- colnames(velo)
cellnames = lapply(strsplit(cellID_velo,':'),'[',2) %>% unlist
cellnames = gsub('x','-1',cellnames)
batch = lapply(strsplit(cellID_velo,':'),'[',1) %>% unlist %>% as.factor()
batch <-  factor(batch, levels = c('S13WT_out','S19WT_out','S31WC_out','S34KO_out','S35KO_out','S7KO_out'),
                 labels = c('_4','_5','_6','_1','_2','_3'))
cellID = paste(cellnames,batch,sep = '')
velo$cellnames <- cellID

velo_sub <- velo[,velo$cellnames %in% colnames(sce)]

velo_sub$Condition <- as.character(sce@meta.data[velo_sub$cellnames,]$Condition) 
velo_sub$celltype <- as.character(sce@meta.data[velo_sub$cellnames,]$celltype_fine)
velo_sub$batch <- as.character(sce@meta.data[velo_sub$cellnames,]$orig.ident)
velo_sub$group <- paste(velo_sub$celltype,velo_sub$batch,sep="_")

### calculate spliced/unspliced proportion
layers_keys = c("spliced", "unspliced", "ambiguous")
counts_layers = sapply(layers_keys,function(key){colSums(velo_sub@assays[[key]]@counts)})
counts_total = rowSums(counts_layers)
counts_total = ifelse(counts_total == 0,counts_total+1,counts_total)
counts_layers = counts_layers/counts_total

data = data.frame(counts_layers)
data$celltype = as.character(velo_sub$celltype)
data$Condition = as.character(velo_sub$Condition)
data$group = paste(data$celltype,data$Condition,sep='_')

df = rbind(aggregate(data$spliced,by=list(celltype=data$group),mean),
           aggregate(data$unspliced,by=list(celltype=data$group),mean),
           aggregate(data$ambiguous,by=list(celltype=data$group),mean))
colnames(df) = c("group" ,"percent")
df$type <- rep(layers_keys,each=length(unique(data$group)))
############################## differentiation potency #############################
library(CytoTRACE)
CytoTRACE_score <- CytoTRACE(mat = as.matrix(Stromal@assays$RNA@counts),enableFast = F,ncores = 5)
Stromal$CytoTRACE_score <- CytoTRACE_score$CytoTRACE[colnames(Stromal)]
################################ monocle2 ######################################
library(monocle)
st = subset(sce,celltype_fine %in% c('Peri','Myo','Fibro'))

expr_matrix <- st@assays$RNA@counts %>% as.matrix() %>% as('sparseMatrix')
p_data <- st@meta.data
f_data <- data.frame(gene_short_name = rownames(st@assays$RNA),row.names = rownames(st@assays$RNA))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds) %>% estimateDispersions()

cds <- detectGenes(cds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
cds <- cds[expressed_genes,]
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds,color_by="celltype_fine", size=2,show_backbone=TRUE)+
plot_cell_trajectory(cds,color_by="Pseudotime", size=2,show_backbone=TRUE)

### find psudotiom gene
Time_diff <- differentialGeneTest(cds, cores = 20,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- filter(Time_diff,qval < 0.05)
Time_genes <- top_n(Time_diff, n = 100, -qval) %>% pull(gene_short_name) %>% as.character()
#################### CellChat all #########################
library(CellChat)
library(ggalluvial)

future::plan("multicore", workers = 20) 
options(future.globals.maxSize = 8000 * 1024^2)

DefaultAssay(sce) <- 'RNA'

cellchatList <- list()
for (i in c('WT','KO')) {
  Object <- subset(sce,(Condition == i))
  cellchat <- createCellChat(object = Object,meta=Object@meta.data,group.by = 'celltype_fine')
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)

  cellchat <- identifyOverExpressedGenes(cellchat) 
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj =PPI.mouse)

  cellchat <- computeCommunProb(cellchat, raw.use = T, population.size=TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  cellchatList[[i]] <- cellchat
}

cellchat <- mergeCellChat(cellchatList, add.names = names(cellchatList))

