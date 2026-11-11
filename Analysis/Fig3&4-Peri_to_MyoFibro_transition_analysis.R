library(ggplot2)
library(ggpubr)
library(Seurat)
library(monocle)
library(UCell)
load("01_rawdata/Internal/Stromal.RData")

DimPlot(Stromal, group.by = "orig.ident")
DimPlot(Stromal, group.by = "seurat_clusters")
DimPlot(Stromal, group.by = "celltype", label = TRUE)

## ---------- Marker expression check ----------
KidneyMarkers <- list(
  PT = c("Lrp2", "Slc27a2", "Pvalb", "Pgam2", "Uqcrb"),
  DCT = c("Slc12a3", "Wnk1"),
  DLH = c("Bst1", "Epha7", "Tmem207"),
  ALH = c("Umod", "Slc12a1", "Car2", "Slc25a4"),
  `CD-IC` = c("Atp6v0d2", "Atp6v1g3", "Fxyd4", "Apela"),
  `CD-PC` = c("Aqp2", "Hsd11b2", "Magi2", "Cdkn1c"),
  Podo = c("Mafb", "Nphs2", "Adgrl4"),
  Endo = c("Plvap", "Emcn", "Tek", "Rgs5"),
  Stromal = c("Col1a1", "Acta2", "Pdgfrb"),
  UE = c("Psca", "Upk3a", "Ly6d"),
  RBC = c("Hbb-bs", "Hba-a1", "Hba-a2")
)

DotPlot(Stromal, features = unique(unname(unlist(KidneyMarkers)))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

markers_for_dotplot <- list(
  Peri = c("Pdgfrb", "Rgs5", "Notch3", "Mylk"),
  Myo = c("Acta2", "Tagln", "Col1a2", "Fn1"),
  Fibro = c("Col1a1", "Col3a1", "Dcn", "Lum")
)
DotPlot(Stromal, features = unique(unname(unlist(markers_for_dotplot)))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

FeaturePlot(
  Stromal,
  features = c("Lrp2", "Slc27a2", "Col1a1", "Acta2"),
  order = TRUE,
  pt.size = 0.5
) |> gs(pic_supp, "FeaturePlot_doublet_check", 8, 8)

stromal_filter <- Stromal[, Stromal$celltype != "Unknow"]
# saveRDS(stromal_filter, "02_processed/scRNA/in_house_stromal_filter.rds")

## ---------- Plot ----------
strcolors = c("Peri" = "#fbdd9a", "Myo" = "#e19d99", "Fibro" = "#99badd")
p <- plot_scanpy_umap(stromal_filter, group_by = "celltype", label = F, custom_colors = strcolors)
legend_only <- cowplot::get_legend(p)
gs(p, pic_fig3, "Fig3B_Stromal_subtypes", width = 5.5, height = 6)
gs(legend_only, pic_fig3, "Fig3B_Stromal_subtypes_legend", width = 5.5, height = 6)

p <- DimPlot(stromal_filter, group.by = "orig.ident", alpha = 1, raster = F) +
  scale_color_manual(values = sample_colors) +
  theme(
    axis.text = element_blank(),  
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14)
  )+
  labs(x="", y="", title = "")
gs(p, pic_supp, "SFig5E_Stromal_batch", width = 4, height = 4)

# deg_stromal <- FindAllMarkers(stromal_filter, group.by = "celltype", only.pos = T,
#                               min.diff.pct = 0.25, min.pct = 0.25, logfc.threshold = 0.25)
# top_stromal <- deg_stromal %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 8)
p <- FacetDotPlot(stromal_filter, markers_for_dotplot, ylab = "Cell Subtype")
gs(p, pic_supp, "SFig4_Stromal_markers", width = 5.5, height = 3.2)

strcolors = c("Peri" = "#fbdd9a", "Myo" = "#e19d99", "Fibro" = "#99badd")
p <- stromal_filter@meta.data %>%
  ggplot(aes(x = reorder(celltype,CytoTRACE_score),y = CytoTRACE_score ,fill = celltype)) +
  geom_violin(trim = F, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  scale_fill_manual(values = strcolors)+
  theme(text = element_text(size = 14),
        axis.text = element_text(color = "black"),
        plot.caption = element_text(size = 8),
        legend.position = "none",
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        panel.background = element_blank()) +
  geom_signif(comparisons = list(c('Myo','Fibro'),c('Peri','Myo'),c('Peri','Fibro')),
              test = "t.test",
              step_increase=0.1, # test.args = c("greater"), 
              map_signif_level = c(
                "*" = 0.05,
                "**" = 0.01,
                "***" = 0.001,
                "****" = 0.0001
              )
              )+
  ylab('CytoTRACE score') +
  xlab("")
gs(p, pic_fig3, "Fig4A_Stromal_stemness", width = 5, height = 5)


## ---------- Differential expression ----------
stromal_filter$celltype_Condition <- paste(stromal_filter$celltype, stromal$Condition)
Idents(stromal_filter) <- "celltype_Condition"

unique(stromal_filter$celltype_Condition)
table(stromal_filter$celltype_Condition)

Peri_deg <- FindMarkers(stromal_filter, ident.1 = "Peri KO", ident.2 = "Peri WT", min.pct = 0.25)
Peri_deg$gene <- rownames(Peri_deg)
p1 <- plot_volcano(Peri_deg, log2fc_col = "avg_log2FC",   
                   pval_col = "p_val", label_col = "gene",
                   upid = "KO", downid = "WT") + theme(legend.position = "none") +
  labs(title = "Peri")

Myo_deg <- FindMarkers(stromal_filter, ident.1 = "Myo KO", ident.2 = "Myo WT", min.pct = 0.25)
Myo_deg$gene <- rownames(Myo_deg)
p2 <- plot_volcano(Myo_deg, log2fc_col = "avg_log2FC", 
                   pval_col = "p_val", label_col = "gene",
                   upid = "KO", downid = "WT") + theme(legend.position = "none") + 
  labs(title = "Myo")

Fibro_deg <- FindMarkers(stromal_filter, ident.1 = "Fibro KO", ident.2 = "Fibro WT", min.pct = 0.25)
Fibro_deg$gene <- rownames(Fibro_deg)
p3 <- plot_volcano(Fibro_deg, log2fc_col = "avg_log2FC", 
                   pval_col = "p_val", label_col = "gene",
                   upid = "KO", downid = "WT") + labs(title = "Fibro")

p <- p1 + p2 + p3
gs(p, pic_fig3,  "Fig3K_stromal_subtype_DEG", width = 13.5, height = 5)

## ---------- ECM score ----------
ECMSigs <- c('Thbs2', 'Col4a4', 'Tnr', 'Col6a6', 'Vwf', 'Gp6', 'Col1a2', 'Itga3', 'Itga9', 'Col4a5', 'Gp5', 'Lamc1', 'Itga6', 'Itga1', 'Itgb5', 'Itgb8', 'Sv2c', 'Thbs3', 'Hmmr', 'Cd44', 'Lama3', 'Lamb1', 'Lamb3', 'Col6a4', 'Lama4', 'Itgb4', 'Lamc3', 'Itgb7', 'Col1a1', 'Fn1', 'Vtn', 'Sv2a', 'Col9a2', 'Itgb3', 'Gp1bb', 'Col9a3', 'Hspg2', 'Npnt', 'Gp1ba', 'Col9a1', 'Itga5', 'Itgav', 'Lama2', 'Col2a1', 'Spp1', 'Sdc4', 'Tnxb', 'Itga10', 'Col4a3', 'Sdc1', 'Lama5', 'Sv2b', 'Col6a3', 'Reln', 'Itga7', 'Itga2', 'Lamc2', 'Col6a2', 'Ibsp', 'Col4a6', 'Comp', 'Col6a1', 'Tnn', 'Itga11', 'Itgb6', 'Lamb2', 'Chad', 'Thbs4', 'Gp9', 'Itga8', 'Dag1', 'Itga2b', 'Agrn', 'Col4a2', 'Col4a1', 'Cd36', 'Lama1', 'Col6a5', 'Itgb1', 'Itga4', 'Thbs1', 'Tnc', 'Cd47')
ECMSigs <- ECMSigs[ECMSigs %in% rownames(sce)]
ECMSigs <- list('ECM' = ECMSigs)
stromal_filter <- AddModuleScore_UCell(stromal_filter, assay='RNA',slot = 'data',features = ECMSigs, ncores = 40, maxRank=2000)
p <- VlnPlot(stromal_filter, group.by = "celltype_Condition", features = "ECM_UCell")

dput(unique(df$ident))
colors <- c("Fibro KO" = "#8e91b3", "Fibro WT" = "#e27979", 
            "Myo KO" = "#8e91b3", "Myo WT" = "#e27979", 
            "Peri KO" = "#8e91b3", "Peri WT" = "#e27979")

p <- ggplot(df, aes(x = ident, y = ECM_UCell, fill = ident)) +
  stat_boxplot(geom="errorbar", width=0.45, size=0.6) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 1) +
  geom_jitter(width = 0.2, size = 1.5) + 
  scale_fill_manual(values = colors) +
  geom_signif(comparisons = list(c("Fibro KO", "Fibro WT"),
                                 c("Myo KO", "Myo WT"),
                                 c("Peri KO", "Peri WT")),
              test = "t.test",
              step_increase = 0,
              map_signif_level = c(
                "*" = 0.05,
                "**" = 0.01,
                "***" = 0.001,
                "****" = 0.0001
              ))
p

## ---------- Functional enrichment ----------
Peri_gene <- Peri_deg %>% filter(p_val < 0.05 & avg_log2FC > 1) %>%  pull(gene)
Peri_fun <- perform_enrichment_analysis(Peri_gene)
p1 <- plot_GO_bubble(
  df = Peri_fun[[1]], 
  topN = 10, font_length = 30, 
  col_vec = c("white", "#e27979", "#d81d1d"),
  title = "Peri"
)

Myo_gene <- Myo_deg %>% filter(p_val < 0.05 & avg_log2FC > 1) %>%  pull(gene)
Myo_fun <- perform_enrichment_analysis(Myo_gene)
p2 <- plot_GO_bubble(
  df = Myo_fun[[1]], 
  topN = 10, font_length = 30, 
  col_vec = c("white", "#e27979", "#d81d1d"),
  title = "Myo"
)

Fibro_gene <- Fibro_deg %>% filter(p_val < 0.05 & avg_log2FC >1) %>% pull(gene) 
Fibro_fun <- perform_enrichment_analysis(Fibro_gene)
p3 <- plot_GO_bubble(
  df = Fibro_fun[[1]], 
  topN = 10, font_length = 30, 
  col_vec = c("white", "#e27979", "#d81d1d"),
  title = "Fibro"
)

p <- p1 + p2 + p3
gs(p, pic_supp, "FigS5M_stromal_upregulated_enrichment", width = 16, height = 4)
## ---------- Monocle2 ----------
DefaultAssay(stromal_filter) <- "RNA"
stromal_monocle <- RunMonocle(stromal_filter)

plot_cell_trajectory(stromal_monocle, color_by = "Pesudotime")+
  plot_cell_trajectory(stromal_monocle, color_by = "celltype") +
  plot_cell_trajectory(stromal_monocle, color_by = "seurat_clusters")

p1 <- plot_cell_trajectory(stromal_monocle, color_by = "celltype", cell_size = 1)  +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
  scale_color_manual(values = strcolors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank(),
    legend.key.size = unit(0.4, "cm"),         
    legend.text = element_text(size = 8),      
    legend.title = element_text(size = 9),      
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),                    
    axis.ticks.length = unit(0.25, "cm")                           
  )+
  labs(color = "Cell Subtype",
       title = "")
p1
gs(p1, pic_fig3, "Fig3B_Stromal_trajectory", width = 5, height = 3.8)

p1 <- plot_cell_trajectory(stromal_monocle, color_by = "Condition", cell_size = 0.3)  +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
  scale_color_manual(values =  c("#d81d1e","#1e3a7e")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank(),
    legend.key.size = unit(0.4, "cm"),        
    legend.text = element_text(size = 8),      
    legend.title = element_text(size = 9),      
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),                   
    axis.ticks.length = unit(0.25, "cm")                          
  )+
  labs(color = "Condition",
       title = "")
gs(p1, pic_fig3, "Fig4B_Stromal_trajectory_condition", width = 5, height = 3.8)

plot_cell_trajectory(stromal_monocle, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters)

# df <- stromal_monocle@phenoData@data
# p <- ggplot(filter(df,celltype == 'Peri'), aes(Pseudotime, colour = Condition, fill=Condition)) +
#   geom_density(bw=0.5,size=1,alpha = 0.5)+
#   ggtitle('Pericyte')+ 
#   theme(text = element_text(size = 8),
#         plot.title = element_text(hjust = 0.5))+
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme_classic()+
#   scale_fill_manual(values = c("#1e3a7e", "#d81d1e"))+
#   scale_color_manual(values = c("#1e3a7e", "#d81d1e"))+
#   theme(legend.position = c(0.8,0.3),legend.direction = 'horizontal')
# gs(p, pic_fig3, "Fig3 inhouse_Peri", width = 7, height = 3.8)

p1 <- plot_cell_trajectory(stromal_monocle, color_by = "Pseudotime", cell_size = 1)  +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
  scale_color_gradientn(colors = c("#1e3a7e", "white", "#d81d1e")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank(),
    legend.key.size = unit(0.4, "cm"),         
    legend.text = element_text(size = 8),      
    legend.title = element_text(size = 9),     
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),                   
    axis.ticks.length = unit(0.25, "cm")                           
  )+
  labs(color = "Pseudotime",
       title = "")
p1
gs(p1, pic_fig3, "Fig4C_Stromal_trajectory_pseudotime", width = 5, height = 3.8)

Time_diff <- differentialGeneTest(stromal_monocle, cores = 20, fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- filter(Time_diff,qval < 0.05)
saveRDS(Time_diff, "03_analysis/Trajectory/inhouse_stromal_time_diff.rds")

df <- data.frame(
  Pseudotime = stromal_monocle$Pseudotime,
  ECM_UCell = stromal_monocle$ECM_UCell,
  celltype = stromal_monocle$celltype,
  state = stromal_monocle$State
)
df <- filter(df, state %in% c(1,3))
p <- ggplot(df, aes(Pseudotime, ECM_UCell, color = celltype)) +
  geom_point(alpha = 1, size = 0.85) + 
  scale_color_manual(values = strcolors)  +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.25, "cm")
  ) +
  labs(color = "Cell Subtype")
gs(p, pic_fig3, "Fig4D ECMscore (inhouse)", 5, 3.8)

Time_genes <- top_n(Time_diff, n = 100, -qval) %>% pull(gene_short_name) %>% as.character()
fibrosis_genes <- c(
  "Col3a1", "Col5a2", "Fn1", "Fbn1", "Timp1", "Bgn", "Dcn",
  "Olfml3", "Col15a1", "Tnc", "Spp1", "Serpine1",
  "Ctgf", "Col1a2", "Fbln2", "Mfap5", "Emp1", "Mgp", "Tgfb1",
  "Cx3cl1", "Cdh11", "Loxl1", "Aebp1", "Sparc",
  "Col1a1",  "Cxcl14", "F2r", "Ltbp2", "Fstl1",
  "Col8a1", "Thbs1",  "S100a6"
)

p <- plot_pseudotime_heatmap(stromal_monocle[fibrosis_genes,],return_heatmap = T,
                        num_clusters = 1,
                        cluster_rows = T,
                        hmcols = colorRampPalette(c("#1e3a7e", "white", "#d81d1e"))(200),
                        cores = 1,
                        show_rownames = T)
gs(p, pic_fig3, "Fig4E inhouse_tarjGenes", height = 6.2, width = 4)


saveRDS(stromal_monocle, "03_analysis/Trajectory/in_house_stromal_monocle.rds")


