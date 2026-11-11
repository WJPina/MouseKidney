## ---------- Packages ----------
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtext)
library(patchwork)
library(reshape2)
library(cowplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(UCell)
library(mclust)
library(ComplexHeatmap)
library(circlize)
library(cellchat)
library(ggpubr)
library(rstatix)
library(ggbreak)

set.seed(12315) 

sample_colors <- c("KO1" = '#D13808', "KO2" = '#FF6836', "KO3"  = '#FFC5AC',
                   "WT1"  = '#227BA2', "WT2" = '#22BAED', "WT3" = '#ABDFFF')


## ---------- Custom functions ----------
source(here("R", "sc_yuanplot.r"))
source(here("R", "yuanplot.r"))

## ---------- Paths ----------
## Figures
fig_dir  <- here("figures")
fig_fig3 <- here("figures", "Fig3")
fig_fig4 <- here("figures", "Fig4")
fig_supp <- here("figures", "Supplementary")
fig_data <- here("figures", "datacollection")

## Data
data_dir <- here("data", "processed", "scRNA")
raw_dir  <- here("data", "raw")

## Create directories if not exist
dir.create(fig_fig3, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_fig4, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_data, recursive = TRUE, showWarnings = FALSE)

## ---------- Plot saving function ----------
gs <- function(plot, path, name, width = 5, height = 5) {
  ggsave(file.path(path, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 600)
  ggsave(file.path(path, paste0(name, ".pdf")),
         plot, width = width, height = height)
}

## ---------- Load data ----------
## Raw (unfiltered) data
seu_raw <- readRDS(here("data", "raw", "Internal", "seu_notfilter.rds"))

## Filtered & integrated data
sce <- readRDS(here("data", "processed", "scRNA", "inhouse_sce_filter.rds"))

## Additional stromal object (used downstream)
load(here("data", "raw", "Internal", "Stromal.RData"))

sample_levels <- c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3")

seu_raw$orig.ident <- factor(seu_raw$orig.ident, levels = sample_levels)
sce$orig.ident     <- factor(sce$orig.ident, levels = sample_levels)

## ---------- Plot ----------
p <- VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0, ncol = 3, cols = sample_colors, group.by = "orig.ident")
gs(p, fig_supp, "SFig5A_QC_before_filtering", width = 9, height = 3.8)

p <- table(seu_raw$orig.ident) |>
  reshape2::melt() |>
  dplyr::rename(
    Sample = Var1, Counts = value
  ) |>
  ggplot(aes(x = Counts, y = Sample, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(color = "black")) +
  geom_text(aes(label = Counts), position = "stack", vjust = 0.5, hjust = 0, size = 3.5)
gs(p, fig_supp, "SFig5B_cell_number_before_filtering", width = 6, height = 3.8)

p <- VlnPlot(sce,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, 
             ncol = 3, cols = sample_colors, group.by = "orig.ident")
gs(p, fig_supp, "SFig4C_QC_after_filtering", width = 9, height = 3.8)

p <- table(sce$orig.ident) |>
  reshape2::melt() |>
  dplyr::rename(
    Sample = Var1, Counts = value
  ) |>
  ggplot(aes(x = Counts, y = Sample, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_text(aes(label = Counts), position = "stack", vjust = 0.5, hjust = 0, size = 3.5) +
  xlim(0, 9000)
gs(p, fig_supp, "SFig4D_cell_number_after_filtering", width = 6, height = 3.8)

p <- DimPlot(sce, group.by = "orig.ident", alpha = 0.7, raster = FALSE) +
  scale_color_manual(values = sample_colors) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + labs(x = "", y = "", title = "")
gs(p, fig_supp, "SFig4E_batch_corrected_UMAP", width = 6, height = 6)

p6 <- plot_scanpy_umap(sce, group_by = "celltype", label = F)
gs(p6, pic_fig3, "Fig3B_UMAP", width = 8, height = 8.5)
legend_only <- cowplot::get_legend(p6)
gs(legend_only, pic_fig3, "Fig3B_UMAP_legend", width = 8, height = 8.5)

KidneyMarkers <- list(PT = c("Slc27a2", "Lrp2"), DCT = c("Pgam2", "Pvalb", "Wnk1"), DLH = c("Epha7", "Bst1"), ALH = c("Tmem207", "Slc12a1", "Umod"), CD_IC = c("Uqcrb", "Slc25a4", "Car2", "Atp6v1g3", "Atp6v0d2"), CD_PC = c("Aqp3", "Apela", "Fxyd4", "Hsd11b2", "Aqp2"), Podo = c("Cdkn1c", "Magi2", "Nphs2", "Mafb"), Endo = c("Tek", "Adgrl4", "Emcn", "Plvap"), Stromal = c("Rgs5", "Pdgfrb", "Acta2", "Col1a1"), UE = c("Ly6d", "Upk3a", "Psca"), RBC = c("Hba-a2", "Hba-a1", "Hbb-bs") )
p <- FacetDotPlot(sce = sce, marker_list = KidneyMarkers)
gs(p, pic_fig3, "Fig3C_Dotplot", width = 12, height = 4.5)

genes <- c("Cdkn1a", "Trp53", "Tgfb1", "Fn1", "Acta2", "Vim")
p <- lapply(genes, function(x)plot_gene_bar(sce, x))
pp <- wrap_plots(p, ncol = 6)
gs(pp, pic_fig3, "Fig3D diffGene", width = 10, height = 3.5)

## ---------- Differential expression ----------
agg_expr <- AggregateExpression(sce, group.by = "orig.ident")[[1]] %>% as.data.frame()
agg_expr <- agg_expr[
  apply(agg_expr, 1, function(x){
    if (sum(x > 0) == 0) return(FALSE)
    x1 <- as.numeric(x[1:3])
    x2 <- as.numeric(x[4:6])
    if (sd(x1) == 0 && sd(x2) == 0) return(FALSE)
    return(TRUE)
  }),
]

#  t-test
res_list <- apply(agg_expr, 1, function(x)t.test(x[4:6], x[1:3]))
p_val <- sapply(res_list,'[[', "p.value")
p_adj <- p.adjust(p_val, method = "BH")
agg_expr$gene <- rownames(agg_expr)
FC <- apply(agg_expr, 1, function(x)mean(x[4:6])/mean(x[1:3]))

deg <- data.frame(genes = names(FC),
                  p_val =p_val, p_adj = p_adj, FC = FC)
saveRDS(deg, "03_analysis/Enrichment/Fig3_all_deg.rds")

## ---------- Functional enrichment ----------
up_genes <- deg_filter$genes[deg_filter$FC > 2 & deg_filter$p_val < 0.05]
down_genes <- deg_filter$genes[deg_filter$FC < 0.5 & deg_filter$p_val < 0.05]
up_fun <- perform_enrichment_analysis(up_genes)
down_fun <- perform_enrichment_analysis(down_genes)

p <- plot_GO_bubble(
  df = up_fun[[1]], 
  topN = 10,
  col_vec = c("white", "#e27979", "#d81d1d")
)
gs(p, fig_fig3, "Fig3E_GO_upregulated", width = 7.5, height = 5.2)
p <- plot_GO_bubble(
  df = down_fun[[1]], 
  topN = 10,
  col_vec = c("white", "#8e91b3", "#1e3a7e")
)
gs(p, fig_supp, "SFig_GO_downregulated", width = 7.5, height = 5.2)

## ---------- Senescence score ----------
AgeSigs = list('Senescence' = AgeSigs)
sce <- AddModuleScore_UCell(sce,
                            assay = 'RNA',slot = 'data',
                            features = AgeSigs, ncores = 40, maxRank = 2000)

## ---------- ECM score ----------
ECMSigs <- c('Thbs2', 'Col4a4', 'Tnr', 'Col6a6', 'Vwf', 'Gp6', 'Col1a2', 'Itga3', 'Itga9', 'Col4a5', 'Gp5', 'Lamc1', 'Itga6', 'Itga1', 'Itgb5', 'Itgb8', 'Sv2c', 'Thbs3', 'Hmmr', 'Cd44', 'Lama3', 'Lamb1', 'Lamb3', 'Col6a4', 'Lama4', 'Itgb4', 'Lamc3', 'Itgb7', 'Col1a1', 'Fn1', 'Vtn', 'Sv2a', 'Col9a2', 'Itgb3', 'Gp1bb', 'Col9a3', 'Hspg2', 'Npnt', 'Gp1ba', 'Col9a1', 'Itga5', 'Itgav', 'Lama2', 'Col2a1', 'Spp1', 'Sdc4', 'Tnxb', 'Itga10', 'Col4a3', 'Sdc1', 'Lama5', 'Sv2b', 'Col6a3', 'Reln', 'Itga7', 'Itga2', 'Lamc2', 'Col6a2', 'Ibsp', 'Col4a6', 'Comp', 'Col6a1', 'Tnn', 'Itga11', 'Itgb6', 'Lamb2', 'Chad', 'Thbs4', 'Gp9', 'Itga8', 'Dag1', 'Itga2b', 'Agrn', 'Col4a2', 'Col4a1', 'Cd36', 'Lama1', 'Col6a5', 'Itgb1', 'Itga4', 'Thbs1', 'Tnc', 'Cd47')
ECMSigs <- ECMSigs[ECMSigs %in% rownames(sce)]
ECMSigs <- list('ECM' = ECMSigs)
sce <- AddModuleScore_UCell(sce, assay='RNA',slot = 'data',features = ECMSigs, ncores = 40, maxRank=2000)

p <- FeaturePlot(
  object = sce,
  features = "ECM_UCell",
  pt.size = 0.3,
  alpha = 1,
  max.cutoff = 0.225, min.cutoff = 0.05,
  order = T
) +
  scale_color_gradientn(colours = c("#83a29a", "#edb760", "#e6807e")) +
  labs(title = "ECM score", x = "UMAP_1", y = "UMAP_2") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_rect(color = "black", fill = NA))
gs(p, fig_fig3, "Fig3G_ECM_score_UMAP", width = 6.5, height = 5)

legend_only <- get_legend(p)
gs(legend_only, fig_fig3, "Fig3G_ECM_score_UMAP_legend")

p <- plot_score_violin(sce,
                       score_col = "ECM_UCell",
                       group_col = "Condition",
                       ylab_text = "ECM Score")
gs(p, fig_fig3, "Fig3F_ECM_score_violin", width = 3.5, height = 5.8)
p <- plot_score_violin(sce,
                       score_col = "Senescence_UCell",
                       group_col = "Condition",
                       ylab_text = "Senescence Score")
gs(p, fig_fig3, "SFig5L_Senescence_score_violin", width = 3.5, height = 5.8)

## ---------- Gaussian mixture model (ECM score) ----------
mycols = c('low'='#84a59d','middle'='#f6bd60','high'='#f28482')

df_box <- data.frame(Celltype = sce$celltype_fine,
                     Condition = sce$Condition,
                     Score = sce$ECM_UCell)

gaussian = Mclust(df_box$Score)
df_box$Group =  factor(ifelse(gaussian$classification %in%c(1,2,3,4),'low',ifelse(gaussian$classification %in% c(5,6,7),'middle','high')),levels = c('low','middle','high'))

p1 <- ggplot(df_box,aes(x= Score,fill = Group))+
  geom_histogram()+
  scale_fill_manual(values = mycols)+
  xlab("ECM Score")+
  ylab('Cell Number')+
  theme_bw()+
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 15),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    text = element_text(size = 20),
    legend.position = 'none'
  )+
  annotate("text", x = 0.07, y = 6000,label = "low",colour=mycols['low'])+
  annotate("text", x = 0.13, y = 2500,label = "middle",colour=mycols['middle'])+
  annotate("text", x = 0.2, y = 500,label = "high",colour=mycols['high'])

sce$ECM_group = df_box[colnames(sce),]$Group
# DimPlot(sce,group.by = 'ECM_group',cols = mycols, split.by = 'Condition')

# Cell-type composition
plot_data <- sce@meta.data %>%
  dplyr::select(celltype_fine, ECM_group) %>%
  group_by(celltype_fine, ECM_group) %>%
  summarise(Count = n(), .groups = 'drop_last') %>%
  mutate(Percent = Count / sum(Count)) %>%
  ungroup()
plot_data <- plot_data[! as.character(plot_data$celltype_fine) %in% c("RBC", "Unknow"),]

ordered_celltypes <- plot_data %>%
  filter(ECM_group == "high") %>%
  arrange(desc(Percent)) %>%
  pull(celltype_fine)

plot_data$celltype_fine <- factor(plot_data$celltype_fine, levels = ordered_celltypes)

p2 <- ggplot(plot_data, aes(x = celltype_fine, y = Percent, fill = ECM_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = mycols) +
  theme_bw() +
  labs(fill = "ECM Group") +
  ylab('Percent') +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    text = element_text(size = 20))

gs(p1, fig_fig3, "Fig3H_GMM_ECM_classification", width = 7, height = 4)
gs(p2, fig_fig3, "Fig3I_ECM_group_by_celltype", width = 8, height = 4)

# VlnPlot(sce, features = "ECM_UCell", pt.size = 0)



## ---------- Cell type proportion ----------
# 1. all 
freq_mat <- as.matrix(table(sce$orig.ident, sce$celltype))
freq_mat <- freq_mat / rowSums(freq_mat)
freq <- reshape2::melt(freq_mat)
colnames(freq) <- c("Batch", "Cluster", "Percent")
freq$Condition <- factor(substr(freq$Batch, 1, 2),
                         levels = c("WT", "KO"))
cluster_order <- freq %>% 
  group_by(Cluster) %>% 
  summarise(mean_percent = mean(Percent)) %>% 
  arrange(desc(mean_percent)) %>% 
  pull(Cluster)
freq$Cluster <- factor(freq$Cluster, levels = cluster_order)

pwc <- freq %>%
  group_by(Cluster) %>%
  t_test(Percent ~ Condition) %>%
  add_significance("p")
pwc <- pwc %>% add_xy_position(x = "Cluster")
pwc$y.position <- pwc$y.position - 0.08
pwc

df_summary <- freq %>%
  group_by(Cluster, Condition) %>%
  summarise(
    mean = mean(Percent),
    sd   = sd(Percent),
    se   = sd(Percent)/sqrt(n()),
    .groups = "drop"
  )

p <- ggplot() +
  geom_col(
    data = df_summary, aes(x = Cluster, y = mean, fill = Condition),
    position = position_dodge(width = 0.6), width = 0.5, alpha = 0.9) +
  geom_errorbar(data = df_summary, aes(x = Cluster, y = mean, ymin = mean - sd, 
                                       ymax = mean + sd, group = Condition),
    width = 0.3, position = position_dodge(width = 0.6), linewidth = 0.6) +
  scale_fill_manual(values = c("KO"="#e68f8b", "WT"="#9fa0be")) +
  geom_point(
    data = freq,aes(x = Cluster, y = Percent, color = Condition),
    position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.12),
    size = 1,alpha = 1) +
  scale_color_manual(values =  c("KO" = "#e7161a", "WT" = "#273b85")) +
  ggpubr::stat_pvalue_manual(
    pwc, hide.ns = TRUE, tip.length = 0, step.increase = 0, label = "p.signif") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),
    legend.key.size = unit(3,'mm'),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = rel(1)),
    panel.background = element_blank()
  ) +
  labs(y = "Percent") +
  scale_y_continuous(limits = c(0, 0.9)) + 
  scale_y_break(breaks = c(0.3,0.6))
gs(p, pic_fig3, "Fig3J_cellratio", width = 7, height = 4)

# 2. stromal 
stromal_filter <- readRDS("02_processed/scRNA/in_house_stromal_filter.rds")
freq_mat = as.matrix(table(stromal_filter$orig.ident,as.character(stromal_filter$celltype)))
freq_mat = freq_mat/rowSums(freq_mat)
freq <- reshape2::melt(freq_mat)
colnames(freq) <- c("Batch", "Cluster", "Percent")
freq$Condition <- factor(substr(freq$Batch, 1, 2),
                         levels = c("WT", "KO"))
cluster_order <- freq %>% 
  group_by(Cluster) %>% 
  summarise(mean_percent = mean(Percent)) %>% 
  arrange(desc(mean_percent)) %>% 
  pull(Cluster)

freq$Cluster <- factor(freq$Cluster, levels = cluster_order)

pwc <- freq %>%
  group_by(Cluster) %>%
  t_test(Percent ~ Condition) %>%
  add_significance("p")
pwc <- pwc %>% add_xy_position(x = "Cluster")
pwc$y.position <- pwc$y.position - 0.08
pwc

df_summary <- freq %>%
  group_by(Cluster, Condition) %>%
  summarise(
    mean = mean(Percent),
    sd   = sd(Percent),
    se   = sd(Percent)/sqrt(n()),
    .groups = "drop"
  )
df_summary$low <- df_summary$mean - df_summary$sd
df_summary$low  <- ifelse(df_summary$low > 0 , df_summary$low , 0)

p <-  ggplot() +
  geom_col(
    data = df_summary, aes(x = Cluster, y = mean, fill = Condition),
    position = position_dodge(width = 0.6), width = 0.5, alpha = 0.9) +
  geom_errorbar(data = df_summary, aes(x = Cluster, y = mean, ymin = mean - sd, 
                                       ymax = mean + sd, group = Condition),
                width = 0.3, position = position_dodge(width = 0.6), linewidth = 0.6) +
  scale_fill_manual(values = c("KO"="#e68f8b", "WT"="#9fa0be")) +
  geom_point(
    data = freq,aes(x = Cluster, y = Percent, color = Condition),
    position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.12),
    size = 1,alpha = 1) +
  scale_color_manual(values =  c("KO" = "#e7161a", "WT" = "#273b85")) +
  ggpubr::stat_pvalue_manual(
    pwc, hide.ns = TRUE, tip.length = 0, step.increase = 0, label = "p.signif") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),
    legend.key.size = unit(3,'mm'),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = rel(1)),
    panel.background = element_blank()
  ) +
  labs(y = "Percent")
gs(p, pic_fig3, "Fig3J_stromal_cellratio", width = 4.2, height = 4)
## ---------- Cell-cell communication (CellChat) ----------
ko <- sce[,sce$Condition == "KO"]
wt <- sce[,sce$Condition == "WT"]
ko <- ko[, ko$celltype_fine != 'Unknow']
wt <- wt[, wt$celltype_fine != 'Unknow']

ko_cellchat <- runCellChat(ko, group.by = "celltype_fine",
                           min.cells = 10, species = "mouse")
saveRDS(ko_cellchat, here("analysis", "CellChat", "ko_cellchat.rds"))

wt_cellchat <- runCellChat(wt, group.by = "celltype_fine",
                           min.cells = 10, species = "mouse")
saveRDS(wt_cellchat, here("analysis", "CellChat", "wt_cellchat.rds"))


WT  <- cc@net[[1]]$weight
KO  <- cc@net[[2]]$weight
diff_mat <- KO - WT
diff_mat[is.na(diff_mat)] <- 0
my_colors <- colorRamp2(
  c(min(diff_mat), 0, max(diff_mat)),  
  c("#1f3a7e", "#ffffff","#d81d1e") 
)
row_sum <- rowSums(diff_mat)
col_sum <- colSums(diff_mat)
row_order <- names(sort(row_sum, decreasing = T))
column_order <-  names(sort(col_sum, decreasing = T))

ha_top <- HeatmapAnnotation(
  Relative_values = anno_barplot(col_sum, border = FALSE,
                        gp = gpar(fill = "#e44b3b"),
                        height = unit(2, "cm"))
)
ha_right <- rowAnnotation(
  Relative_values = anno_barplot(row_sum, border = FALSE,
                        gp = gpar(fill = "#e44b3b"),
                        width = unit(2, "cm"))
)

pdf(here("figures", "Fig4", "Fig4L_KO_vs_WT_CellChat_heatmap.pdf"),
    width = 6.5, height = 4.8)
ht <- Heatmap(
  diff_mat, col = my_colors, name = "Relative values",
  row_order = row_order, column_order = column_order,
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = ha_top, right_annotation = ha_right,
  
  row_names_side = "left", row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = 9), row_title_rot = 90,
  
  column_title = "Differential interaction strength (KO vs WT)", column_title_gp = gpar(fontsize = 10)
)
draw(ht)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = "Myo",
                                            color.use = c("grey10","#1e3a7e","#d81d1e"))

pdf(here("figures", "Fig4", "Fig4M_Myo_signaling_changes.pdf"), width = 7, height = 5)
gg1
dev.off()

gg2 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = "Peri",
                                            color.use = c("grey10","#1e3a7e","#d81d1e"))

pdf(here("figures", "Fig4", "Fig4N_Peri_signaling_changes.pdf"), width = 7, height = 5)
gg2
dev.off()

gg3 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = "Fibro",
                                            color.use = c("grey10","#1e3a7e","#d81d1e"))

pdf(here("figures", "Fig4", "Fig4_Fibro_signaling_changes.pdf"), width = 7, height = 5)
gg3
dev.off()

