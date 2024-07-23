library(ggplot2)
library(ggrepel)
library(dplyr)
library(magrittr)
library(tibble)
library(Seurat)
library(rstatix)
library(ggpubr)
library(scCustomize)
library(introdataviz)
library(reshape2)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

mypalette <- read.csv("./colors.csv",header = T)

setwd('./Figures/')

################################## QC ########################################
### plot raw quality
RawMerge <- merge(Raw[[1]], Raw[c(2:6)])
RawMerge$percent.mt <- PercentageFeatureSet(RawMerge, pattern = "^mt-")
RawMerge$orig.ident <- factor(RawMerge$orig.ident,levels = c("WT1","WT2","WT3","KO1","KO2","KO3"))

mycols = c('#D13808','#FF6836','#FFC5AC','#227BA2','#22BAED','#ABDFFF')
names(mycols) <- c("WT1","WT2","WT3","KO1","KO2","KO3")

png("RawQC.png",width = 2000,height = 800,res = 300)
VlnPlot(RawMerge, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), pt.size = 0,
        ncol = 3,cols = mycols,group.by = 'orig.ident')
dev.off()

### counts
png( "Batch_counts.png", width = 2000,height = 1800,res = 300)
table(sce$orig.ident) %>% 
  melt() %>%
  rename("Frequncey"= "value") %>%
  rename("Batch"="Var1") %>%
  ggplot(aes(x = Frequncey,y = Batch,fill = Batch)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = mycols) + 
  theme_classic()+
  theme(text = element_text(size = 18)) +
  labs(fill = "")+
  geom_text(aes(label = Frequncey),position="stack", vjust=0.5,hjust=1,col="black",size=5)
dev.off()

### Condition
png("All_Condition.png",width = 2500,height = 2000,res = 300)
DimPlot_scCustom(sce, figure_plot = TRUE,colors_use = c("#007E99","#FF745A"),group.by = 'Condition',pt.size = 0.5)
dev.off()


### Cluster
png("All_Cluster.png",width = 2500,height = 2000,res = 300)
DimPlot_scCustom(sce, figure_plot = TRUE,group.by = 'seurat_clusters',label = T,colors_use = mypalette$palette3,
                 label.box = T,repel = T,pt.size = 0.5)&theme(legend.position = 'none')
dev.off()

### Cluster
png( "Cluster_percent.png", width = 2500,height = 2000,res = 300)
table(sce$orig.ident,sce$seurat_clusters) %>% 
  melt() %>%
  rename("Frequncey"= "value") %>%
  rename("Batch"="Var1") %>%
  rename("cluster"="Var2") %>% 
  ggplot(aes(y = Frequncey,x = Batch,fill = reorder(cluster,Frequncey))) +
  geom_bar(stat = "identity",position = "fill")+
  scale_fill_manual(values = mypalette$palette3) + 
  theme_classic()+
  theme(text = element_text(size = 18)) +
  labs(fill = "")
dev.off()
################################## Anno ########################################
Idents(sce) <- 'celltype'
mycols <- c("#C71585","#EAA86C","#F7DDD4","#FFE4B5","#827141","#992C73","#7673AE","#EAA944","#74B5CF","#E18772","#CA9600")
names(mycols) <- levels(sce)
### umap
png("All_umap.png",width = 2500,height = 2000,res = 300)
DimPlot_scCustom(sce,figure_plot = TRUE,colors_use = mycols,label.box = T,repel = T,pt.size = 0.5)&theme(legend.position = 'none')
dev.off()

png("All_AnnoDot.png",width = 1400,height = 350,res = 200)
table(sce$seurat_clusters) %>% 
  as.data.frame() %>% 
  mutate(Celltype = c('PT', 'PT', 'PT', 'DLH', 'PT', 'PT', 'PT', 'PT', 'Endo', 'ALH', 'Endo', 'PT', 'PT', 'DCT', 'PT', 
                      'Stromal', 'Endo', 'CD-IC', 'CD-PC', 'PT', 'PT', 'Podo', 'UE', 'Stromal', 'DLH', 'RBC')) %>%
  rename('Var1'='Cluster') %>%
  rename('Freq'='Cell number') %>%
  ggplot(aes(x=Cluster,y=1,size=`Cell number`,color = Celltype))+
  geom_point()+
  scale_size(range = c(2,10))+
  scale_color_manual(values = mycols)+theme_classic()+
  theme(axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.line.x = element_blank(),legend.position = 'bottom',text = element_text(size=15))+
  guides(color = F)
dev.off()

png("All_umap_condition.png",width = 2000,height = 1000,res = 300)
DimPlot(sce, cols = mycols,pt.size = 0.5,split.by = 'Condition')
dev.off()

### marker
# png("All_markerDot.png",width = 2000,height = 2500,res = 200)
  pdf("All_markerDot.pdf",width = 12,height = 6)
  DotPlot(sce,features = signatures,scale = T,assay = "RNA",group.by = 'celltype') +
    scale_color_gradient2(low = "blue",high = "red",midpoint = 0) + 
    theme(axis.text.x = element_text(size = 15,angle = 90,vjust = 0.5,hjust = 1),
          axis.title.x = element_blank(),
          axis.line = element_line(colour = "black",linewidth = rel(1)))
  dev.off()
################################## PT  ########################################
### cluster
png("ST_umap_cluster.png",width = 1200,height = 1200,res = 300)
DimPlot(Stromal, label.box = T,repel = T,pt.size = 0.5,group.by = 'seurat_clusters',
  label = T,cols=c("#984EA3","#E41A1C","#377EB8","#4DAF4A"))&
  theme(legend.position = 'none')&
  ggtitle('Stromal clusters')
dev.off()

### ref
mycols <- mypalette$palette9[1:4]
names(mycols) <- c("S1",'S2','S3','S3T2')

Idents(ref_use) <- 'celltype'
p <- DimPlot_scCustom(ref_use,  colors_use = mycols,figure_plot = TRUE,label.box = T,repel = T,pt.size = 0.25)&theme(legend.position = 'none')
p[[1]] <- p[[1]]+ggtitle('PT reference')
png("PT_umap_ref.png",width = 1000,height = 800,res = 300)
p
dev.off()

### transfer
Idents(PT) <- 'celltype_pre'
p <- DimPlot_scCustom(PT, figure_plot = TRUE,colors_use = mycols,label.box = T,repel = T,pt.size = 0.25)&theme(legend.position = 'none')
p[[1]] <- p[[1]]+ggtitle('PT transfer labels')
png("PT_umap_trans.png",width = 1000,height = 1000,res = 300)
p
dev.off()

### correlation
png("PT_umap_refcor.png",width = 1200,height = 1200,res = 250)
corrplot(cor_mat,method="square",col = rev(COL2('RdBu', 10)),addCoef.col = 'black',col.lim = c(-0.2,0.7),
         tl.col = 'black', tl.srt = 0,tl.offset = 0.8,
         is.corr = F,order = 'AOE',addrect = 3,rect.col = 'red', rect.lwd = 2)
dev.off()

### umap cluster annotation
Idents(PT) <- 'celltype'
mycols <- mypalette$palette9[1:4]
names(mycols) <- levels(PT$celltype)

p <- DimPlot_scCustom(PT, figure_plot = TRUE,colors_use = mycols,label.box = T,repel = T,pt.size = 0.25)&
  theme(legend.position = 'none')
p[[1]] <- p[[1]]+ggtitle('PT subtypes')
png("PT_umap.png",width = 1000,height = 1000,res = 300)
p
dev.off()
################################## Stromal  ###################################
### cluster
png("ST_umap_cluster.png",width = 1200,height = 1200,res = 300)
DimPlot(Stromal, label.box = T,repel = T,pt.size = 0.5,group.by = 'seurat_clusters',
  label = T,cols=c("#984EA3","#E41A1C","#377EB8","#4DAF4A"))&
  theme(legend.position = 'none')&
  ggtitle('Stromal clusters')
dev.off()

### marker
Idents(Stromal) <- 'seurat_clusters'
pList <- list()
for(i in names(signatures_ST)){
  pList[[i]] <- Stacked_VlnPlot(Stromal, features = signatures_ST[[i]], x_lab_rotate = TRUE,colors_use = c("#984EA3","#E41A1C","#377EB8","#4DAF4A"))
}
png("ST_marker.png",width = 1600,height = 800,res = 250)
ggarrange(pList[[1]],pList[[2]],pList[[3]],ncol = 3,nrow = 1,common.legend = F,widths = c(4,4,4))
dev.off()

### umap
mycols <- mypalette$palette3[1:4]
names(mycols) <- levels(Stromal$celltype)
Idents(Stromal) <- 'celltype'

p <- DimPlot_scCustom(Stromal, figure_plot = TRUE,colors_use = mycols,label.box = T,repel = T,pt.size = 1)&theme(legend.position = 'none')
p[[1]] <- p[[1]]+ggtitle('Stromal subtypes')
png("ST_umap.png",width = 1000,height = 1000,res = 300)
p
dev.off()

### heatmap
DEGs <- FindAllMarkers(Stromal,only.pos = T,return.thresh = 0.05)
Topmarker <- DEGs[order(DEGs$cluster,DEGs$avg_log2FC,decreasing = T),] %>% group_by(cluster) %>% do(head(., n = 20))

group <- Stromal@active.ident
my_order <- sort(table(Stromal@active.ident),decreasing = T)
group <- group[order(factor(group,levels = names(my_order)))]
group <- factor(group,levels = names(my_order))

Topmarker <- Topmarker[order(factor(Topmarker$cluster,levels=names(my_order))),] 
Stromal <- ScaleData(Stromal,features = Topmarker$gene)
mat <- as.matrix(Stromal@assays$RNA@scale.data[Topmarker$gene,names(group)])

col <- colorRamp2(c(-1,0,1), c("lightgrey","white", "red"), transparency = 0.3, space = "LAB")
top_anno <- HeatmapAnnotation(df = data.frame(CellType = group),
                              col = list(CellType = mycols),
                              show_legend = T,
                              annotation_legend_param = list(title = "",
                                                             labels_gp = gpar(fontsize = 12), 
                                                             title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                             nrow = 1))
h <- Heatmap(mat,
             show_column_names = F,
             show_row_names = T,
             col = col,
             cluster_rows = F,
             cluster_columns = F,
             top_annotation = top_anno,
             row_names_gp = gpar(fontsize = 12),
             heatmap_legend_param = list(title = "Expression\nLevel",
                                         at = c(-1,0,1),
                                         labels = c("-1", "0","1"),
                                         labels_gp = gpar(fontsize = 15),
                                         title_gp = gpar(fontsize = 15, fontface = "bold"),
                                         legend_width = unit(30, "mm")))
png("ST_heatmap.png",width = 2000,height = 2500,res = 200)
draw(h,heatmap_legend_side = "right",annotation_legend_side = "top",merge_legend = F,align_heatmap_legend = "heatmap_center")
dev.off()
########################## celltype fine ###################################
Idents(sce) <- 'celltype_fine'
mycols = c("#FF4500", "#1AAF8B", "#406C85", "#F6BD16", 
           "#EAA86C","#F7DDD4","#FFE4B5","#827141","#992C73","#7673AE","#EAA944",
           "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
           "#E18772","#CA9600")
names(mycols) = levels(sce$celltype_fine)
### celltype condition
png("All_umap_condition.png",width = 2500,height = 1200,res = 300)
DimPlot(sce, cols = mycols,pt.size = 0.5,split.by = 'Condition')
dev.off()

png( "Celltype_percent.png", width = 1500,height = 1800,res = 300)
table(sce$Condition,sce$celltype_fine) %>% 
  melt(value.name = 'Frequncey') %>%
  rename("Var1"="Condition") %>%
  rename("Var2"="Celltype") %>% 
  ggplot(aes(y = Frequncey,x = Condition,fill = reorder(Celltype,Frequncey))) +
  geom_bar(stat = "identity",position = "fill")+
  scale_fill_manual(values = mycols) + 
  theme_classic()+
  theme(text = element_text(size = 18)) +
  ylab("Percent")+
  guides(fill = guide_legend(title = 'Celltype'))
dev.off()

########################## senescence phenotype #######################
### Vdr
sig = t.test(subset(sce['Vdr',],Condition=='WT')@assays$RNA@data %>% as.numeric,
             subset(sce['Vdr',],Condition=='KO')@assays$RNA@data %>% as.numeric,
             alternative = c("greater"))
sig

png("All_Vdr.png",width = 1500,height = 2000,res = 300)
(FeaturePlot(subset(sce,Condition=='WT'),'Vdr',order = T,pt.size = 0.5)+ggtitle('WT'))/
  (FeaturePlot(subset(sce,Condition=='KO'),'Vdr',order = T,pt.size = 0.5)+ggtitle('KO',)+labs(caption = 'One-tail t.test p-value < 2.2e-16'))
dev.off()

### senescence/ECM markers
genes <- c("Cdkn1a","Trp53","Tgfb1","Fn1","Acta2","Vim")

png("All_seneGenes.png",width = 3000,height = 2500,res = 300)
data.frame(t(as.matrix(sce@assays$RNA@data[genes,])),Condition = sce$Condition) %>%
  melt(value.name = 'Expression Level',variable.name ='Gene',id.name='Condition') %>% 
  ggplot(aes(x = Condition,y=`Expression Level`,fill=Condition)) +
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  geom_violin()+
  geom_jitter(size=0.1,width = 0.25)+
  facet_wrap(~Gene,scales = 'free_y')+
  theme_classic2()+
  theme(strip.background = element_blank(),text = element_text(size=25),
        legend.position = 'none',axis.title.x = element_blank())+
  geom_signif(comparisons = list(c("KO","WT")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("greater")) 
dev.off()

### Senescence score
df_box <- data.frame(Celltype = sce$celltype,
                     Condition = sce$Condition,
                     Score = sce$Senescence_UCell)

pwc <- df_box %>% 
        t_test(Score  ~ Condition) %>% 
        adjust_pvalue(method = "bonferroni") %>% 
        add_significance("p.adj")
pwc <- pwc %>% add_xy_position(x = "Condition")
pwc

png("Senescence_UCell.png",width = 900,height = 1000,res = 300)
df_box %>%
  ggplot() +
  geom_boxplot(aes(x = Condition,y = Score ,fill = Condition))+
  ggpubr::stat_pvalue_manual(pwc,hide.ns = T) +
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  theme(text = element_text(size = 15),strip.background = element_blank(),
        plot.caption = element_text(size = 8),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank()) +
  ylab('Senescence Score')
dev.off()
#################### plot differential genes #########################
### DEG volcano
res <- DEGs_Con
res$threshold = ifelse(res$p_val_adj <= 0.05, ifelse(res$avg_log2FC > 0,"Up","Down"),"Background")

png("DEG_Condition.png",width = 1500,height = 2000,pointsize = 24,res = 150)
ggplot(data = res, aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold)) +
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c("Up" = "#FF745A","Down"="#007E99","Background" = "#ABABAB"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="grey50",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey50",lwd=0.6) +
  labs(title = "WT vs KO")+
  ylab('-Log10 (P-value)')+
  xlab('LogFC')+
  theme_bw()+
  xlim(c(-4,4))+
  scale_x_continuous(breaks = seq(-4,4,1))+
  theme(
    plot.title = element_text(size=15,hjust = 0),
    axis.title.x = element_text(size=15,),
    axis.title.y = element_text(size=15,),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank())
dev.off()

### DEGs heatmap
for(i in c('down','up')){
  if (i == 'down'){
    genes <- filter(res,avg_log2FC < -1 & p_val_adj == 0)
  }
  else {genes <- filter(res,avg_log2FC >1 & p_val_adj == 0)}
  
  Object <- sce
  genes = rownames(genes)
  Object <- ScaleData(Object,features = genes)
  mat = Object@assays$RNA@scale.data[genes,]
  mat = t(mat)
  matmean = aggregate(mat,by = list(Condition = Object$Condition),mean)
  matmean  <- matmean %>% column_to_rownames('Condition') %>% t()
  col <- circlize::colorRamp2(c(-0.5,0,0.5), c("#007E99","white", "#FF745A"), space = "LAB")
  h <- Heatmap(matmean,
          show_column_names = T,
          show_row_names = T,
          rect_gp = gpar(col = "lightgrey", lty = 1, lwd = 0.5),
          col = col,
          cluster_rows = T,
          cluster_columns = F,
          row_names_gp = gpar(fontsize = 12),
          heatmap_legend_param = list(title = "mean\nexpression"))
  if(i == 'down'){
    hei  = 1800
    wid = 900
  }
  else{
    hei=2400
    wid = 800}
  png(paste("Mean_heatmap_",i,".png",sep=''),width = wid,height = hei,res = 250)
  print(h)
  dev.off()
}
#################### plot Enrichment between Condition #########################
p1<-
df_up[order(df_up$Count,decreasing = T),][1:10,] %>%
  ggplot(aes(reorder(Description, Count),Count,fill = qvalue))+
  geom_col()+
  theme_bw()+
  scale_fill_gradient(low = "#D72E3E",high = "#FEB1A3")+
  theme(legend.position = 'bottom',
        legend.key.size = unit(12,'mm'),
        legend.key.height = unit(3,'mm'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(color="black",size=15),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  coord_flip()+
  ylim(0, 20)+
  geom_text(data = df_up[order(df_up$Count,decreasing = T),][1:10,] ,aes(label = reorder(Description, Count),y= 0.01),size = 3,hjust = 0)+
  guides(size = FALSE)+
  ggtitle('Up-regulated genes (KO vs WT)')+ theme(plot.title = element_text(hjust = 0))+
  xlab('GO pathways')+ylab('Gene counts')

p2<-
df_down[order(df_down$Count,decreasing = T),][1:10,] %>%
  ggplot(aes(reorder(Description, Count),Count,fill = qvalue))+
    geom_col()+
    theme_bw()+
    scale_fill_gradient(low = "#2E5A87",high = "#B0C0D2")+
    theme(legend.position = 'bottom',
          legend.key.size = unit(12,'mm'),
          legend.key.height = unit(3,'mm'),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text = element_text(color="black",size=15),
          axis.line.x = element_line(color='black'),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())+
    coord_flip()+
    ylim(0, 20)+
    geom_text(data = df_down[order(df_down$Count,decreasing = T),][1:10,],aes(label = reorder(Description, Count),y= 0.01),size = 3,hjust = 0)+
    guides(size = FALSE)+
    ggtitle('Down-regulated genes (KO vs WT)')+ theme(plot.title = element_text(hjust = 0))+
    xlab('GO pathways')+ylab('Gene counts')

png("GO.png",width = 2500,height = 1500,res = 350)
ggarrange(p1,p2,ncol = 2)
dev.off()

#################### ECM phenotype #########################
library(mclust)
mycols = c('low'='#84a59d','middle'='#f6bd60','high'='#f28482')

df_box <- data.frame(Celltype = sce$celltype_fine,
                     Condition = sce$Condition,
                     Score = sce$ECM_UCell)
set.seed(1234)
gaussian = Mclust(df_box$Score)
df_box$Group =  factor(ifelse(gaussian$classification %in%c(1,2,3,4),'low',ifelse(gaussian$classification %in% c(5,6,7),'middle','high')),levels = c('low','middle','high'))

png("ECM_UCell.png",width = 1200,height = 1000,res = 300)
ggplot(df_box,aes(x= Score,fill = Group))+
  geom_histogram()+
  scale_fill_manual(values = mycols)+
  xlab("ECM Score")+
  ylab('Cell Number')+
  theme_bw()+
  theme(axis.text.x =element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 12),legend.position = 'none')+
  annotate("text", x = 0.07, y = 6000,label = "low",colour=mycols['low'])+
  annotate("text", x = 0.13, y = 2500,label = "middle",colour=mycols['middle'])+
  annotate("text", x = 0.2, y = 500,label = "high",colour=mycols['high'])
dev.off()

sce$ECM_group = df_box[colnames(sce),]$Group

png("ECM_UCell_Con.png",width = 2800,height = 1500,res = 300)
DimPlot(sce,group.by = 'ECM_group',cols = mycols,split.by = 'Condition')
dev.off()

Object <- subset(sce,celltype=='Stromal')
freq = as.matrix(table(Object$orig.ident,Object$ECM_group))
freq = freq/rowSums(freq)
freq = data.frame(freq)
colnames(freq) <- c("Batch","Cluster","Percent")
freq$Condition <- factor(substr(freq$Batch,1,2))
freq = filter(freq,Cluster == 'high')
freq$Condition = factor(freq$Condition,levels = c('WT','KO'))

png("ST_ECM_UCell_high_percent.png",width = 1000,height = 1500,res = 300)
freq %>%
ggplot(aes(x=Condition,y=Percent,fill=Condition)) +
  geom_boxplot()+
  geom_jitter()+
  geom_signif(comparisons = list(c("WT","KO")),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("less"))+
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  theme_bw()+
  theme(axis.text.x =element_text(size = 15,angle = 90,vjust = 0.5,hjust = 1),
        legend.position = 'none',
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 20))+
  ylab('Percent of high ECM group')
dev.off()

png("ECM_UCell_percent.png",width = 2200,height = 2000,res = 300)
ggplot(sce@meta.data) +
  geom_bar(aes(x=celltype_fine,fill=ECM_group),position='fill',stat='Count')+
  scale_fill_manual(values = mycols)+
  theme_bw()+
  theme(axis.text.x =element_text(size = 15,angle = 90,vjust = 0.5,hjust = 1),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 20))+
  ylab('Percent')
dev.off()

### Stromal ECM
Stromal$ECM_UCell <- sce@meta.data[colnames(Stromal),'ECM_UCell']

png("ST_ECM_UCell.png",width = 1200,height = 1000,res = 300)
FeaturePlot_scCustom(Stromal,features = 'ECM_UCell',label = T,pt.size = 1.5)+ggtitle('ECM Score')
dev.off()

############################# Cell type proportion ########################
library(gg.gap)

png("proportion.png",width = 2000,height = 1600,res = 300)
p <-
freq %>%
  mutate(Condition = factor(Condition,levels=c('WT','KO'))) %>%
  ggplot()+ 
  geom_bar(aes(x = Cluster,y=Percent,fill = Condition),stat = "identity",width = 0.5,position = "dodge")+ 
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  ggpubr::stat_pvalue_manual(pwc,hide.ns = T,tip.length = 0,step.increase = 0,label =  "p.signif")+
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)))

gg.gap(plot = p,
       segments = c(0.2, 0.7),
       tick_width = 0.1,
       rel_heights = c(0.2, 0, 0.1),
       ylim = c(0,0.9))

dev.off()


png("ST_proportion.png",width = 900,height = 800,res = 400)
freq %>%
  mutate(Condition = factor(Condition,levels=c('WT','KO'))) %>%
  ggplot()+ 
  geom_bar(aes(x = Cluster,y=Percent,fill = Condition),stat = "identity",width = 0.5,position = "dodge")+ 
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  theme_bw()+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 0,hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(3,'mm'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank())+
  ggpubr::stat_pvalue_manual(pwc,hide.ns = T,tip.length = 0,step.increase = 0,label =  "p.signif")
dev.off()

################################### RNA velocity #################################
### RNA spliced proportion
png(paste("Spliced_proportion_all.png",sep = ""),width = 1400,height = 2300,res = 300)
ggplot(df,aes(y=group,x=percent,fill=type))+
  geom_bar(stat = 'identity',position = 'fill')+
  theme_classic()+
  theme(plot.margin = margin(20,0,5,0),
        text = element_text(size = 15),axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.position = c(0.38,1.02),legend.key.size = unit(3,'mm'),legend.key.height = unit(3,'mm'),legend.key.width = unit(3,'mm'),
        legend.text = element_text(size=10),legend.direction = 'horizontal',legend.title = element_text(size=10,vjust = 1))+
  geom_text(aes(label = ifelse(type== 'unspliced', paste0(sprintf("%.0f", percent*100),"%"),"")),position="stack",
            vjust=0.5,hjust=0,col="white",size=5)+
  scale_fill_manual(values = c('spliced'="#0a9396",'ambiguous'="#e9d8a6",'unspliced'="#9b2226"))
dev.off()


png(paste("Spliced_proportion_subset.png",sep = ""),width = 900,height = 1500,res = 300)
filter(df,group %in% c('Peri_KO','Peri_WT','Myo_KO','Myo_WT','Fibro_KO','Fibro_WT')) %>%
ggplot(aes(x=group,y=percent,fill=type))+
  geom_bar(stat = 'identity',position = 'fill')+
  theme_classic()+
  theme(plot.margin = margin(20,0,5,0),
        text = element_text(size = 15),axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.position = c(0.38,1.05),legend.key.size = unit(3,'mm'),legend.key.height = unit(3,'mm'),legend.key.width = unit(3,'mm'),
        legend.text = element_text(size=10),legend.direction = 'horizontal',legend.title = element_text(size=10,vjust = 1))+
  geom_text(aes(label = ifelse(type== 'unspliced', paste0(sprintf("%.0f", percent*100),"%"),"")),position="stack",
            vjust=0,hjust=0.4,col="white",size=4)+
  scale_fill_manual(values = c('spliced'="#0a9396",'ambiguous'="#e9d8a6",'unspliced'="#9b2226"))
dev.off()
################################### differentiation potency #################################
mycols <- mypalette$palette3[1:4]
names(mycols) <- levels(Stromal$celltype)

png("ST_CytoTRACE_score.png",width = 800,height = 1200,res = 300)
Stromal@meta.data %>%
  ggplot(aes(x = reorder(celltype,CytoTRACE_score),y = CytoTRACE_score ,fill = celltype)) +
  geom_boxplot()+
  scale_fill_manual(values = mycols)+
  theme(text = element_text(size = 14),
        plot.caption = element_text(size = 8),
        legend.position = "none",
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank()) +
  geom_signif(comparisons = list(c("Unknow","Fibro"),c('Myo','Fibro'),c('Peri','Myo')),
              test = "t.test",
              step_increase=0.1,
              map_signif_level = T,
              test.args = c("greater"))+
  xlab('celltype')
dev.off()


Stromal$group = paste(Stromal$celltype,Stromal$Condition,sep='_')

p1 <- ggplot(filter(Stromal@meta.data,celltype == 'Peri') )+
  geom_density(aes(x = CytoTRACE_score,color = group,fill = group),alpha = 0.3)+
  theme_classic()+
  scale_color_manual(values = c('Peri_WT'="#007E99",'Peri_KO'="#FF745A"))&
  theme(text = element_text(size = 15),
        plot.caption = element_text(size = 8),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank())
p2 <-  ggplot(filter(Stromal@meta.data,celltype == 'Peri'))+
  geom_density(aes(x = ECM_UCell,color = group,fill = group),alpha = 0.3)+
  theme_classic()&
  scale_color_manual(values = c('Peri_WT'="#007E99",'Peri_KO'="#FF745A"))&
  theme(text = element_text(size = 15),
        plot.caption = element_text(size = 8),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank())
png("Peri_score.png",width = 900,height = 1500,res = 300)
ggarrange(p1,p2,ncol = 1,common.legend = T)
dev.off()
#################################### plot trajectory ###########################
mycols <- mypalette$palette3[1:4]
names(mycols) <- levels(Stromal$celltype)

p1<-plot_cell_trajectory(cds,color_by="Pseudotime", size=4,show_backbone=TRUE)+
  theme(legend.position = 'top',
        legend.key.size = unit(5,'mm'),
        legend.key.height = unit(3,'mm'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10,vjust = 1),
        text = element_text(size = 10),
        axis.line = element_line(colour = "black",linewidth = rel(1)))
p2<-plot_cell_trajectory(cds,color_by="celltype_fine", size=4,show_backbone=TRUE)+
  scale_color_manual(values = mycols)+
  theme(legend.position = 'top',
        legend.key.size = unit(4,'mm'),
        legend.key.height = unit(3,'mm'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10,vjust = 1),
        text = element_text(size = 10),
        axis.line = element_line(colour = "black",linewidth = rel(1)))
png(paste("ST_trajectory.png",sep = ""),width = 1800,height = 900,res = 300)
ggarrange(p1,p2,ncol = 2)
dev.off()

library(ggpubr)
df <- cds@phenoData@data
png(paste("ST_trajectory_density.png",sep = ""),width = 2000,height = 1000,res = 300)
ggplot(filter(df,celltype_fine == 'Peri'), aes(Pseudotime, colour = Condition, fill=Condition)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+
  ggtitle('Pericyte')+ 
  theme(text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  scale_color_manual(values = c("#007E99","#FF745A"))+
  theme(legend.position = c(0.8,0.15),legend.direction = 'horizontal')
dev.off()

png(paste("ST_psudotiom_gene_heatmap.png",sep = ""),width = 1200,height = 2400,res = 300)
plot_pseudotime_heatmap(cds[Time_genes,],return_heatmap = T,
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
dev.off()

#################### cellchat #########################
mycols = c("#007E99","#FF745A")
png("cellchat_all_pathway.png",width = 1500,height = 2700,res = 300)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1, 2),measure = "weight",color.use = mycols)
dev.off()
