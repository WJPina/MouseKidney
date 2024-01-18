library(ggplot2)
library(ggrepel)
library(dplyr)
library(magrittr)
library(tibble)
library(Seurat)
library(rstatix)
library(introdataviz)
library(ggpubr)

mypalette <- read.csv("~/scripts/colors.csv",header = T)
################################## UMAP ########################################
### CT
png("Figures/All_umap_CT1.png",width = 2500,height = 2000,res = 300)
DimPlot(sce, reduction = "umap",label = F,pt.size = 0.7,repel = T,label.box = T)+
  scale_color_manual(values = mypalette$palette2)+
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 25),
        legend.position = "bottom",
        legend.text = element_text(size = 12,face = 'bold'),
        legend.key.size = unit(2,"mm"),
        axis.title = element_blank())+
  geom_text_repel(data = aggregate(sce@reductions$umap@cell.embeddings,by = list(celltype = sce$celltype_L2),mean),
                  aes(label= celltype, x= UMAP_1,y= UMAP_2))+
  guides(color=guide_legend(nrow=1,override.aes = list(size = 8)))
dev.off()

### Condition
png("Figures/All_Condition.png",width = 2500,height = 2000,res = 300)
DimPlot(sce, reduction = "umap",label = F,pt.size = 0.4,group.by = "Condition",cols = c("#007E99","#FF745A"))+
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 25),
        legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank())+
  geom_text_repel(data = aggregate(sce@reductions$umap@cell.embeddings,by = list(celltype = sce$celltype_L2),mean),
                  aes(label= celltype, x= UMAP_1,y= UMAP_2))
dev.off()
################################### mean heatmap ###############################
Object = sce
DefaultAssay(Object) <- "RNA"
Object <- ScaleData(Object,features = markers)
mat = Object@assays$RNA@scale.data[markers,]
mat = t(mat)
matmean = aggregate(mat,by = list(celltype = Object@active.ident),mean)
matmean  <- matmean %>% column_to_rownames('celltype') %>% t()
col <- circlize::colorRamp2(c(-1,0,1), c("lightgrey","white", "red"), transparency = 0.3, space = "LAB")
library(ComplexHeatmap)
h = Heatmap(matmean,show_column_names = T,
            show_row_names = T,
            col = col,
            cluster_rows = F,
            cluster_columns = F,
            row_names_gp = gpar(fontsize = 12),
            heatmap_legend_param = list(title = "mean expression",
                                        title_position = 'lefttop',
                                        legend_direction = "horizontal",
                                        at = c(-1,0,1),
                                        labels = c( "-1", "0","1"),
                                        labels_gp = gpar(fontsize = 15),
                                        title_gp = gpar(fontsize = 15, fontface = "bold"),
                                        legend_width = unit(20, "mm")))
png("Figures/Mean_heatmap.png",width = 1000,height = 2500,res = 300)
draw(h,heatmap_legend_side = "top",merge_legend = F,align_heatmap_legend = "heatmap_left")
dev.off()

################################## markers dot #################################
markers = unique(unlist(signatures[levels(sce$celltype_L2)]))
png("Figures/All_markerDot.png",width = 1300,height = 2000,res = 200)
DotPlot(sce,features = markers,scale = T,assay = "RNA",group.by = 'celltype_L2') + 
  scale_color_gradient2(low = "blue",high = "red",midpoint = 0) + 
  RotatedAxis()+
  coord_flip() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)))
dev.off()

########################### plot Cell type proportion ##########################
### t-test
png("Figures/proportion.png",width = 1000,height = 800,res = 300)
freq %>%
  mutate(Condition = factor(Condition,levels=c('WT','KO'))) %>%
  ggplot()+ 
  geom_bar(aes(x = Cluster,y=Percent,fill = Condition),stat = "identity",width = 0.5,position = "dodge")+ 
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  ggpubr::stat_pvalue_manual(pwc,hide.ns = T,tip.length = 0,step.increase = 0,label =  "p.signif")+
  theme(text = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.position = c(0.9,0.89),
        legend.key.size = unit(3,'mm'),
        legend.title = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)))
dev.off()


png("Figures/proportion_Stroma.png",width = 600,height = 800,res = 300)
freq %>%
  mutate(Condition = factor(Condition,levels=c('WT','KO'))) %>%
  ggplot()+ 
  geom_bar(aes(x = Cluster,y=Percent,fill = Condition),stat = "identity",width = 0.5,position = "dodge")+ 
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 0,hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = c(0.2,0.9),
        legend.key.size = unit(3,'mm'),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank())+
  ggpubr::stat_pvalue_manual(pwc,hide.ns = T,tip.length = 0,step.increase = 0,label =  "p.signif")
dev.off()
################################### plot  score ################################
df_box <- data.frame(celltype = sce$celltype_L1,
                     Condition = sce$Condition,
                     group = sce$celltype_L1,
                     Senescence_UCell = sce$Senescence_UCell,
                     ECM_UCell = sce$ECM_UCell)

df_box$score = df_box$Senescence_UCell

### Senesce
pwc <- df_box %>% 
  t_test(score  ~ Condition) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj")
pwc <- pwc %>% add_xy_position(x = "Condition")
pwc

png("Figures/Senescence_UCell.png",width = 600,height = 800,res = 300)
df_box %>%
  ggplot() +
  geom_boxplot(aes(x = Condition,y = score ,fill = Condition))+
  ggpubr::stat_pvalue_manual(pwc,hide.ns = T) +
  scale_fill_manual(values = c("#007E99","#FF745A"))+
  theme(text = element_text(size = 10),
        plot.caption = element_text(size = 8),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank()) +
  ylab('Senescence_UCell')
dev.off()

### ECM
df_box$score = df_box$ECM_UCell

png("Figures/ECM_UCell.png",width = 1000,height = 800,res = 300)
ggplot(df_box,aes(x= celltype,y= score,fill= Condition))+
  geom_split_violin(trim= T,color="black",scale = "width",show.legend = F)+
  scale_fill_manual(values = c('#007E99','#FF745A'))+
  ylab("ECM_UCell")+
  stat_compare_means(data = df_box,
                     aes(x=celltype,y=score,group=Condition),
                     hide.ns = T,
                     label.y = c(0.17,0.24,0.19,0,0.29,0.19,0.24),
                     label = "p.signif")+
  theme(text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth = rel(1)))
dev.off()

#################### plot Enrichment between Condition #########################
col = ifelse(res$avg_log2FC[1]>0,"Reds","Blues")
png("Figures/Peri.GO_down.png",width = 1500,height = 1400,res = 300)
ggplot(dfPlot,aes(reorder(Description, Count),Count,fill = qvalue))+
  geom_col()+
  theme_bw()+
  scale_fill_distiller(palette=col,direction =-1)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  coord_flip()+
  # scale_y_discrete(expand = expansion(mult = c(2,2)))+
  ylim(0, 21)+
  geom_text(data = dfPlot,aes(label = reorder(Description, Count),y= 0.01),size = 3,hjust = 0)+
  guides(size = FALSE)+
  ggtitle('Peri: KO vs WT down-regulated genes')+ theme(plot.title = element_text(hjust = 0))+
  xlab('GO pathways')+ylab('Enrichment gene counts')
dev.off()

#################################### plot trajectory ###########################
png(paste("Figures/Stroma_trajectory.png",sep = ""),width = 1800,height = 1400,res = 300)
plot_cell_trajectory(cds,color_by="Pseudotime", size=4,show_backbone=TRUE)+
  theme(legend.key.height = unit(3,'mm'),
        legend.position = 'top',
        legend.title = element_text(vjust = 1),
        text = element_text(size = 20),
        axis.line = element_line(colour = "black",linewidth = rel(1)))
dev.off()

library(ggpubr)
df <- cds@phenoData@data
png(paste("Figures/Stroma_trajectory_density_KO.png",sep = ""),width = 1200,height = 1000,res = 300)
ggplot(filter(df,Condition == 'KO'), aes(Pseudotime, colour = celltype_L2, fill=celltype_L2)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+
  ggtitle('KO')+ 
  theme(text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c('#61DDAA','#5AB8DB','#FFC771'))+
  scale_color_manual(values = c('#61DDAA','#5AB8DB','#FFC771'))+
  theme(legend.position = 'bottom')
dev.off()

png(paste("Figures/Stroma_psudotiom_gene_heatmap.png",sep = ""),width = 1200,height = 1600,res = 300)
plot_pseudotime_heatmap(cds[Time_genes,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
dev.off()

png(paste("Figures/Stroma_psudotiom_gene_heatmap.png",sep = ""),width = 1200,height = 1600,res = 300)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-100)),],
                            branch_point = 1,
                            num_clusters = 2,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()
