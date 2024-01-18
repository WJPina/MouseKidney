library(Seurat)
library(DoubletFinder)
library(dplyr)
library(magrittr)
library(ggplot2)
############################## Read raw data ###################################
setwd("~/wangj2/MouseKidney/counts/")
folders=list.files('.')
Raw = lapply(folders,function(folder){CreateSeuratObject(counts = Read10X_h5(folder), project = substr(folder,1,3),min.cells = 3,min.features = 200)})
setwd('/home/wangjing2/wangj2/MouseKidney/RData6')

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
############################# Integrate Data ###################################
### Integrate
features <- SelectIntegrationFeatures(object.list = NormList)
Anchors <- FindIntegrationAnchors(object.list = NormList,anchor.features = features,dims = 1:npcs)
Integrated <- IntegrateData(anchorset = Anchors, dims = 1:npcs) 
DefaultAssay(Integrated) <- "integrated"

### Find variable genes from expression level bins
source('~/wangj2/MouseKidney/MouseKineyPro/FindVarBins.R')
Integrated$project <- 'Integrated'
Integrated <- FindVarBins(Integrated)

npcs = 30
Integrated <- ScaleData(Integrated, verbose = FALSE) %>%
              RunPCA(npcs = npcs, verbose = FALSE) %>%
              RunUMAP(dims = 1:npcs,reduction = 'pca') %>%
              FindNeighbors(dims = 1:npcs, reduction="pca")
ElbowPlot(Integrated,ndims = npcs)
DimPlot(Integrated,reduction = 'umap',group.by = "Condition")
########################## Annotate cell type ##################################
sce <- Integrated_var
DefaultAssay(sce)<-'integrated'
sce$Condition <- factor(sce$Condition,levels = c('WT','KO'))
sce$orig.ident <- factor(sce$orig.ident,levels = c('WT1','WT2','WT3','KO1','KO2','KO3'))
sce <- FindClusters(sce,algorithm = 2,resolution = 0.8)
DimPlot(sce,reduction = 'umap',label = T,repel = T,label.box = T,group.by = 'seurat_clusters')+NoLegend()

DefaultAssay(sce)<-'RNA'
signatures = read.csv("/mnt/data3/wangj2/MouseKidney/MarkerList/combineMarkers.csv",header = T,na.strings = "") %>%
  as.list() %>%
  lapply(function(x){x <- x[!is.na(x) & x %in% rownames(sce@assays$RNA)]})
signatures = signatures[lengths(signatures)!=0]
signatures

DotPlot(sce,features = signatures,scale = T,assay = "RNA")+
  scale_color_gradient2(low = "white",high = "blue",midpoint = 0) +
  RotatedAxis()

celltype <- as.character(sce$integrated_snn_res.0.8)
celltype[celltype %in% 21] <- 'UE'
celltype[celltype %in% 23] <- 'Podo'
celltype[celltype %in% c(10)] <- 'ALH'
celltype[celltype %in% c(6,24)] <- 'DLH'
celltype[celltype %in% c(18)] <- 'DCT'
celltype[celltype %in% c(13)] <- 'CD-PC'
celltype[celltype %in% c(19)] <- 'CD-IC'
celltype[celltype %in% c(12,16,9)] <- 'Endo'
celltype[celltype %in% c(15,22)] <- 'Stromal'
celltype[celltype %in% c(7,8)] <- 'PT.S1' 
celltype[celltype %in% c(0,2,3,17,14)] <- 'PT.S2' 
celltype[celltype %in% c(1,4,5,11,20)] <- 'PT.S3' 

sce$celltype <- factor(celltype,levels = c('UE','Podo','ALH','DLH','DCT','CD-PC','CD-IC','Endo','Stromal','PT.S1','PT.S2','PT.S3'))
Idents(sce) <- 'celltype'
DimPlot(sce,label = T)

############################ Recluster Stromal #################################
library(reticulate)
library(sceasy)
library(Seurat)
use_python("/home/wangjing/miniconda3/envs/scvi-env/bin/python")
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
set.seed(123)

Stromal_raw <- subset(sce,celltype %in% c('Stromal'))
Stromal_raw$Cluster <- droplevels(Stromal_raw$seurat_clusters)
DefaultAssay(Stromal_raw) <- 'RNA'
Stromal_raw <- SCTransform(Stromal_raw,vars.to.regress = c('nCount_RNA','nFeature_RNA','percent.mt'))
Stromal_raw <- CreateSeuratStromal_raw(counts = Stromal_raw@assays$SCT@data,meta.data = Stromal_raw@meta.data)
adata <- convertFormat(Stromal_raw, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)

pcs = 10
scvi$settings$seed = as.integer(123)
scvi$model$SCVI$setup_anndata(adata,batch_key = "orig.ident")
model = scvi$model$SCVI(adata,n_latent=as.integer(pcs))
scvi$model$SCVI$view_anndata_setup(model)
model$train(early_stopping = T,early_stopping_monitor = 'train_loss')

latent = model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(Stromal_raw)
Stromal_raw[["scvi"]] <- CreateDimReducStromal_raw(embeddings = latent, key = "scvi_", assay = DefaultAssay(Stromal_raw))
Stromal_raw <- RunUMAP(Stromal_raw, dims = 1:pcs, reduction = 'scvi')
Stromal_raw <- FindNeighbors(Stromal_raw, dims = 1:pcs, reduction = 'scvi')
Stromal_raw <- FindClusters(Stromal_raw,algorithm = 2,resolution = 0.2)

DimPlot(Stromal_raw,label = T)
FeaturePlot(Stromal_raw,features = c('nCount_SCT'),order = T,label = T)

### annotated subclsuter
VlnPlot(Stromal_raw,features = c('Pdgfrb','Acta2','Rgs5','Dcn','Gpx3','Fxyd2','Lrp2'),stack = T)

DEGs <- FindAllMarkers(Stromal_raw,only.pos = T,return.thresh = 0.05)

celltype <- as.character(Stromal_raw$seurat_clusters)
celltype[celltype == 0] <- 'LowPT'
celltype[celltype == 1] <- 'Peri'
celltype[celltype == 2] <- 'Myo'
Stromal_raw$celltype <- factor(celltype,levels = c('Peri','LowPT','Myo'))
Idents(Stromal_raw) <- 'celltype'
DimPlot(Stromal_raw,label = T)

### reclsuter
Stromal <- subset(Object,celltype %in% c('Peri','Myo'))
DimPlot(Stromal,label = T)

adata <- convertFormat(Stromal, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)

pcs = 10
scvi$model$SCVI$setup_anndata(adata,batch_key = "orig.ident")
model = scvi$model$SCVI(adata,n_latent=as.integer(pcs))
scvi$model$SCVI$view_anndata_setup(model)
model$train(early_stopping = T,early_stopping_monitor = 'train_loss')

latent = model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(Stromal)
Stromal[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(Object))
Stromal <- RunUMAP(Stromal, dims = 1:pcs, reduction = 'scvi')
Stromal <- FindNeighbors(Stromal, dims = 1:pcs, reduction = 'scvi')
Stromal <- FindClusters(Stromal,algorithm = 2,resolution = 0.2)

DimPlot(Stromal,label = T)
FeaturePlot(Stromal,features = c('Pdgfrb','Acta2','Rgs5','Dcn'),order = T,label = T)

celltype <- as.character(Stromal$seurat_clusters)
celltype[celltype == 0] <- 'Peri'
celltype[celltype == 1] <- 'Myo'
Stromal$celltype <- factor(celltype,levels = c('Peri','Myo'))
Idents(Stromal) <- 'celltype'
DimPlot(Stromal,label = T)

Stromal <- ScaleData(Stromal)
DEGs <- FindAllMarkers(Stromal,only.pos = T,return.thresh = 0.05)
TopDEGs <- DEGs[order(DEGs$cluster,DEGs$avg_log2FC,decreasing = T),] %>% group_by(cluster) %>% do(head(., n = 20))
DoHeatmap(Stromal,features = TopDEGs$gene )
############################## Cell type proportion #############################
### t-test
library(ggpubr)
library(rstatix)
# freq = as.matrix(table(sce$orig.ident,sce@active.ident))
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
################################ DEGs ##########################################
Object <- sce
DefaultAssay(Object) <- "RNA"
Object<- ScaleData(Object,features = rownames(Object))

### All
Idents(Object) <- 'celltype'
levels(Object@active.ident)
DEGs <- FindAllMarkers(Object,only.pos = T,return.thresh = 0.05,assay = "RNA",min.diff.pct = 0.2,logfc.threshold = 0.25)

### Condition
Idents(Object) <- "Condition"
DEGs_Con <- FindMarkers(Object,only.pos = F,return.thresh = 0.05,assay = "RNA",ident.1 = "KO")

### Condition between type
Idents(Object) <- "celltype"
clusters = levels(Object@active.ident)
DEGs_ConClus = lapply(clusters, function(x){
  print(x)
  if ( x %in% clusters) {
    DEG_ConClu <- FindMarkers(Object,ident.1 = "KO",group.by = "Condition",subset.ident = x,only.pos = F,return.thresh = 0.05,assay = "RNA")
  }
})
names(DEGs_ConClus) <- clusters
sapply(DEGs_ConClus,dim)

save(DEGs,DEGs_Con,DEGs_ConClus,file = 'DEGs.RData')

library(openxlsx)
write.xlsx(DEGs_ConClus,file = 'DEGs_ConClus.xlsx',rowNames = T)
# ################################ Enrichment ####################################
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
# res <- filter(DEGs_Con,avg_log2FC < -1 & p_val_adj < 0.05)
res = filter(DEGs_ConClus$Peri,avg_log2FC < -0.5 & p_val_adj < 0.05)
dim(res)
genelist <- bitr(rownames(res), fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)
bg=bitr(rownames(sce@assays$RNA),'SYMBOL','ENTREZID','org.Mm.eg.db')

### GO
ego <- enrichGO(gene = genelist, OrgDb = 'org.Mm.eg.db',universe = bg$ENTREZID,ont = "All")
go <- ego@result[,c("Description","qvalue","GeneRatio",'Count')]
go$Ratio <- as.numeric(unlist(lapply(strsplit(go$GeneRatio,"/"),"[",1)))/as.numeric(unlist(lapply(strsplit(go$GeneRatio,"/"),"[",2)))
dim(go)
dfPlot <- go %>% filter(qvalue < 0.05)

write.xlsx( GoEnrich,file = 'GOEnrich.xlsx',rowNames = T)


dfPlot <- dfPlot[c('GO:1901615', 'GO:0006631', 'GO:0044282', 'GO:0016324', 'GO:0045177', 'GO:0046394', 'GO:0016053', 'GO:0015711', 'GO:0005777', 'GO:0042579'),]
head(dfPlot)
# ##################################### score ####################################
library(UCell)
DefaultAssay(sce) <- "RNA"
### Aging score
AgeSigs <- read.csv("~/wangj2/GeneSets/MouseAgingList.csv",header = T)
AgeSigs <- AgeSigs$Symbol
AgeSigs = AgeSigs[AgeSigs %in% rownames(sce)]
AgeSigs = list('Senescence' = AgeSigs)
sce <- AddModuleScore_UCell(sce,assay='RNA',slot = 'data',features=AgeSigs,ncores = 40,maxRank=2000)
# ### ECM
ECMSigs =c('Thbs2', 'Col4a4', 'Tnr', 'Col6a6', 'Vwf', 'Gp6', 'Col1a2', 'Itga3', 'Itga9', 'Col4a5', 'Gp5', 'Lamc1', 'Itga6', 'Itga1', 'Itgb5', 'Itgb8', 'Sv2c', 'Thbs3', 'Hmmr', 'Cd44', 'Lama3', 'Lamb1', 'Lamb3', 'Col6a4', 'Lama4', 'Itgb4', 'Lamc3', 'Itgb7', 'Col1a1', 'Fn1', 'Vtn', 'Sv2a', 'Col9a2', 'Itgb3', 'Gp1bb', 'Col9a3', 'Hspg2', 'Npnt', 'Gp1ba', 'Col9a1', 'Itga5', 'Itgav', 'Lama2', 'Col2a1', 'Spp1', 'Sdc4', 'Tnxb', 'Itga10', 'Col4a3', 'Sdc1', 'Lama5', 'Sv2b', 'Col6a3', 'Reln', 'Itga7', 'Itga2', 'Lamc2', 'Col6a2', 'Ibsp', 'Col4a6', 'Comp', 'Col6a1', 'Tnn', 'Itga11', 'Itgb6', 'Lamb2', 'Chad', 'Thbs4', 'Gp9', 'Itga8', 'Dag1', 'Itga2b', 'Agrn', 'Col4a2', 'Col4a1', 'Cd36', 'Lama1', 'Col6a5', 'Itgb1', 'Itga4', 'Thbs1', 'Tnc', 'Cd47')
ECMSigs = ECMSigs[ECMSigs %in% rownames(sce)]
ECMSigs = list('ECM' = ECMSigs)
sce <- AddModuleScore_UCell(sce,assay='RNA',slot = 'data',features=ECMSigs,ncores = 40,maxRank=2000)
################################ monocle2 ######################################
library(monocle)
Object = subset(sce,celltype_L1 == 'Stroma')
expr_matrix <- Object@assays$RNA@counts %>% as.matrix() %>% as('sparseMatrix')
p_data <- Object@meta.data
f_data <- data.frame(gene_short_name = rownames(Object@assays$RNA),row.names = rownames(Object@assays$RNA))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds) %>% estimateDispersions()

Object <- FindVariableFeatures(Object ,assay = 'RNA',selection.method = "vst", nfeatures = 1000)
ordergene <- VariableFeatures(Object,assay = 'RNA')

cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree',norm_method = 'log')
cds <- orderCells(cds)

plot_cell_trajectory(cds,color_by="celltype", size=2,show_backbone=TRUE)
 
### find psudotiom gene
ordergene = cds@featureData@data$gene_short_name[cds@featureData@data$use_for_ordering]
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 20,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- top_n(Time_diff, n = 50, -qval) %>% pull(gene_short_name) %>% as.character()

BEAM_res <- BEAM(cds, branch_point = 1, cores = 20)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]






CellCycleScoring()
