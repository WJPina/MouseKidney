## ---------- Packages ----------
library(Seurat)
library(ggplot2)
library(Matrix)
library(tidyr)
library(cowplot)
library(dplyr)
library(monocle)
library(here)
library(UCell)
library(CytoTRACE)
library(ggpubr)

set.seed(12135)

## ---------- Custom functions ----------
source(here("R", "sc_yuanplot.r"))

## ---------- Paths ----------
fig_dir  <- here("figures")
fig_fig4 <- here("figures", "Fig4")
fig_supp <- here("figures", "Supplementary")

data_dir <- here("data", "processed", "scRNA")
raw_dir  <- here("data", "raw", "GEO")

dir.create(fig_fig3, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_fig4, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp, recursive = TRUE, showWarnings = FALSE)

## ---------- Plot saving function ----------
gs <- function(plot, path, name, width = 5, height = 5) {
  ggsave(file.path(path, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 600)
  ggsave(file.path(path, paste0(name, ".pdf")),
         plot, width = width, height = height)
}

## ---------- Stromal markers ----------
markers_for_stromal <- list(
  VDR   = c("VDR"),
  Peri  = toupper(c("Pdgfrb", "Rgs5", "Notch3", "Mylk")),
  Myo   = toupper(c("Acta2", "Tagln", "Col1a2", "Postn")),
  Fibro = toupper(c("Col1a1", "Col3a1", "Dcn", "Lum"))
)

## ---------- Kidney cell-type markers (detailed) ----------
kidney_markers <- list(
  
  ## --- Proximal tubule ---
  PT_S1 = c("SLC5A2", "LRP2", "SLC34A1", "AQP1", "SLC5A12"),
  PT_S2 = c("SLC22A8", "SLC13A3", "SLC22A6", "SLCO1A1", "NAT8"),
  PT_S3 = c("SLC7A13", "SLC27A2", "SLC7A12", "CPN1", "ATP11A"),
  PT_injured = c("VCAM1", "HAVCR1", "LCN2", "VIM", "SOX9"),
  
  ## --- Loop of Henle ---
  Thin_DL = c("AQP1", "SLC14A2"),
  Thin_AL = c("CLDN10", "KCNJ1"),
  TAL     = c("UMOD", "SLC12A1", "CLDN16", "FXYD2", "KCNJ1"),
  
  ## --- Distal nephron ---
  DCT = c("SLC12A3", "PVALB", "TRPM6", "CALB1"),
  CNT = c("AQP2", "SLC8A1", "CALB1", "SCNN1G"),
  
  ## --- Collecting duct ---
  CD_PC  = c("AQP2", "SCNN1G", "AVPR2", "AQP3"),
  CD_ICA = c("ATP6V1B1", "SLC4A1", "AE1"),
  CD_ICB = c("SLC26A4", "ATP6V1G3"),
  
  ## --- Glomerulus ---
  Podocyte = c("NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO"),
  PEC = c("CLDN1", "KRT8", "KRT18", "LYPD5"),
  Mesangial = c("PDGFRB", "RGS5", "MYH11", "COL4A1", "COL4A2"),
  Endo_Glomerular = c("PLVAP", "EMCN", "ESM1", "KDR", "FLT1"),
  Endo_Peritubular = c("FLT1", "KDR", "VWF", "PECAM1"),
  
  ## --- Stromal ---
  Fibroblast = c("DCN", "LUM", "COL1A1", "COL3A1", "FBLN1"),
  Pericyte   = c("PDGFRB", "RGS5", "NOTCH3", "KCNJ8", "ABCC9"),
  Myofibroblast = c("ACTA2", "TAGLN", "COL1A1", "COL1A2", "FN1", "POSTN"),
  
  ## --- Immune ---
  T_cells = c("CD3D", "CD3E", "CD4", "CD8A"),
  B_cells = c("MS4A1", "CD79A"),
  Plasma_cells = c("MZB1", "IGHG1"),
  Macrophage = c("LYZ", "C1QA", "C1QB", "MSR1"),
  Dendritic = c("FCER1A", "LAMP3")
)

## ---------- Kidney markers (simplified for DotPlot) ----------
kidney_markers_filter <- list(
  PT   = c("SLC5A2", "LRP2", "SLC22A8", "SLC13A3", "SLC27A2"),
  TAL  = c("UMOD", "SLC12A1", "CLDN16", "KCNJ1"),
  CD_IC = c("ATP6V1B1", "SLC4A1", "AE1", "SLC26A4", "ATP6V1G3"),
  CD_PC = c("AQP2", "SCNN1G", "AVPR2", "AQP3"),
  Podocyte = c("NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO"),
  Endo_ = c("EMCN", "KDR", "FLT1", "CALCRL"),
  Pericyte = c("PDGFRB", "RGS5", "ACTA2", "TAGLN", "COL1A1"),
  T_cells = c("CD3D", "CD3E", "CD2", "CD8A"),
  B_cells = c("MS4A1", "CD79A", "MZB1", "IGHG1"),
  Myeloid = c("LYZ", "C1QA", "C1QB", "FCER1A")
)



## ---------- Dataset 1: GSE211785 ----------
GSE211785 <- readRDS(file.path(data_dir, "GSE211785_seu_obj_scRNA_normal.rds"))

GSE211785$project   <- "GSE211785"
GSE211785$cell_type <- GSE211785$Cluster_Idents
GSE211785$orig.ident <- GSE211785$orig_ident

## ---------- Dataset 2: GSE202109 ----------
GSE202109 <- readRDS(file.path(data_dir, "GSE202109_combined_QC.rds"))
GSE202109 <- JoinLayers(GSE202109)

new.cluster.id <- setNames(rep(NA, 17), as.character(0:16))
new.cluster.id[c("0","1","2","3","4","5","6","7","9","10")] <- "PT"
new.cluster.id["8"]  <- "TAL"
new.cluster.id["11"] <- "Stroma"
new.cluster.id["12"] <- "Lymphocytes"
new.cluster.id["13"] <- "CD_PC"
new.cluster.id["14"] <- "CD_ICA"
new.cluster.id["15"] <- "Myeloid"
new.cluster.id["16"] <- "CD_ICB"

GSE202109$cell_type <- unname(new.cluster.id[as.character(GSE202109$seurat_clusters)])
GSE202109$project   <- "GSE202109"

## ---------- Dataset 3: GSE183276 ----------
counts <- readRDS(file.path(
  raw_dir, "GSE183276",
  "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS"
))

meta <- vroom::vroom(file.path(
  raw_dir, "GSE183276",
  "GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt"
))

GSE183276 <- CreateSeuratObject(counts = counts, project = "GSE183276")
GSE183276 <- AddMetaData(GSE183276, metadata = meta)
colnames(GSE183276@meta.data)[22] <- "group"

GSE183276 <- subset(GSE183276, subset = group == "Normal Reference")
GSE183276$project <- "GSE183276"

GSE183276[["RNA"]] <- split(GSE183276[["RNA"]], f = GSE183276$orig.ident)
GSE183276 <- NormalizeData(GSE183276) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

GSE183276 <- IntegrateLayers(
  GSE183276,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)

GSE183276 <- FindNeighbors(GSE183276, reduction = "harmony") %>%
  FindClusters() %>%
  RunUMAP(reduction = "harmony", dims = 1:30)


## ---------- Dataset 4: GSE140989 ----------
GSE140989 <- readRDS(file.path(data_dir, "GSE140989_ann.rds"))
GSE140989$project <- "GSE140989"

new.cluster.id <- setNames(rep(NA, 24), as.character(0:23))
new.cluster.id[c("0","1")] <- "PT"
new.cluster.id["2"] <- "CD_PC"
new.cluster.id[c("3","9","23")] <- "T cells"
new.cluster.id[c("4","10","16")] <- "TAL"
new.cluster.id["5"] <- "CD_ICA"
new.cluster.id["19"] <- "CD_ICB"
new.cluster.id[c("6","15")] <- "Endo_Glomerular"
new.cluster.id["18"] <- "Endo_Arterial"
new.cluster.id["21"] <- "Mesangial"
new.cluster.id[c("7","13")] <- "CNT"
new.cluster.id["8"] <- "PT_injured"
new.cluster.id["11"] <- "Stroma"
new.cluster.id[c("12","17")] <- "Myeloid"
new.cluster.id["14"] <- "DCT"
new.cluster.id["20"] <- "Podo"
new.cluster.id["22"] <- "B cells"

GSE140989$cell_type <- unname(new.cluster.id[as.character(GSE140989$seurat_clusters)])

## ---------- Merge public datasets ----------
seu_obj <- merge(GSE202109, list(GSE140989, GSE183276_nor, GSE211785_sc_nor))
seu_obj <- JoinLayers(seu_obj)
# Batch inspection
seu_obj <- NormalizeData(seu_obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(seu_obj, ndims = 50)
seu_obj <- RunUMAP(seu_obj, dims = 1:30)
DimPlot(seu_obj, reduction = "umap", group.by = "cell_type", label = T, repel = T) + NoLegend()
DimPlot(seu_obj, reduction = "umap", group.by = "project", label = T, repel = T) + NoLegend()

# Harmony integration
seu_obj[["RNA"]] <- split(seu_obj[["RNA"]], f = seu_obj$project)
seu_obj <- NormalizeData(seu_obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(seu_obj, ndims = 50)
seu_obj <- RunUMAP(seu_obj, dims = 1:30)
seu_obj <- IntegrateLayers(seu_obj, method = HarmonyIntegration,
                           orig.reduction = "pca", new.reduction = "harmony")
seu_obj <- FindNeighbors(seu_obj, reduction = "harmony") %>% FindClusters() %>% 
  RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "uamp.harmony")

DimPlot(seu_obj, reduction = "uamp.harmony", group.by = "project")
DimPlot(seu_obj, reduction = "uamp.harmony", group.by = "orig.ident") + NoLegend()
DimPlot(seu_obj, reduction = "uamp.harmony", group.by = "cell_type")
DimPlot(seu_obj, reduction = "uamp.harmony", group.by = "seurat_clusters", label = T) + NoLegend()

# saveRDS(seu_obj, file.path(data_dir, "all_scRNA.rds"))

## ---------- cell-type annotation  ----------
seu_obj <- readRDS(file.path(data_dir, "all_scRNA.rds"))
seu_obj <- FindClusters(seu_obj)
DimPlot(seu_obj, reduction = "uamp.harmony", group.by = "seurat_clusters", label = T) + NoLegend()
FacetDotPlot(seu_obj, kidney_markers, group_by = "seurat_clusters")

new.cluster.id <- setNames(rep(NA, 28), as.character(0:27))
new.cluster.id[c("0", "1", "2", "3", "10","22", "27")] <- "PT"
new.cluster.id[c("8", "9")]  <- "TAL"
new.cluster.id[c("8", "9")]  <- "TAL"
new.cluster.id[c("12", "25")]  <- "Stromal"
new.cluster.id[c("7", "14", "19","21")]  <- "CD-PC"
new.cluster.id[c("13", "17")]  <- "CD-IC"
new.cluster.id[c("23")]  <- "Podo"
new.cluster.id[c("4", "5", "24")]  <- "T cells"
new.cluster.id[c("15")]  <- "B cells"
new.cluster.id[c("11", "16")]  <- "Myeloid"
new.cluster.id[c("6", "18", "20","26")]  <- "Endo"

seu_obj$cell_type_refine <- unname(new.cluster.id[as.character(seu_obj$seurat_clusters)])
seu_obj$cell_type_refine <- factor(seu_obj$cell_type_refine,
                                   levels = c("PT", "TAL", "CD-IC", "CD-PC", "Podo", "Endo", "Stromal", "T cells", "B cells", "Myeloid"))

## ---------- Visualization (Supplementary Fig 6 A-C) ----------
colors <- c("PT" = "#a0d6cc", "Stromal" = "#e07988", "Endo" = "#2b6068", "Podo" = "#f0d6d1", "TAL" = "#d75da6",
            "CD-IC" = "#89c590", "CD-PC" = "#579ab8", "T cells" = "#e7756a", "B cells" = "#f2c387", "Myeloid" = "#9b83b2")
p <- plot_scanpy_umap(seu_obj, group_by = "cell_type_refine", label = F, reduction = "uamp.harmony",custom_colors = colors)
gs(p, fig_dir, "SFig6B_public_data_UMAP", 8, 9)

legend_only <- get_legend(p)
gs(legend_only, fig_dir, "SFig6B_public_data_UMAP_legend", 8, 9)

p <- FacetDotPlot(seu_obj, kidney_markers_filter, group_by = "cell_type_refine", ylab = "Cell Type", 
             colors = c("white", "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#43A2CA"))
gs(p, fig_dir, "SFig6C_marker_dotplot", 13, 4.5)

p <- DimPlot(
  seu_obj, reduction = "uamp.harmony", group.by =  "project", 
  label = F, repel = TRUE, split.by = "project",
  label.size = label_size,  raster = T
) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    legend.position = "top", 
    legend.text = element_text(size = 10), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(x="", y="", title = "") 

gs(p, fig_dir, "SFig6A_batch_mixing", 20, 6)
legend_only <- get_legend(p)
gs(legend_only, fig_dir, "SFig6A_batch_mixing_legend", 8, 9)

## ---------- Stromal subtype annotation ----------
Fibro <- seu_obj[, seu_obj$seurat_clusters %in% c("12","25")] 
Fibro <- NormalizeData(Fibro) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
ElbowPlot(Fibro, ndims = 50) 
Fibro <- RunUMAP(Fibro, dims = 1:30) 
DimPlot(Fibro, group.by = "project") 

Fibro[["RNA"]] <- as(Fibro[["RNA"]], Class = "Assay5") 
Fibro[["RNA"]] <- split(Fibro[["RNA"]], f = Fibro$project) 
Fibro <- NormalizeData(Fibro) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
Fibro <- IntegrateLayers(Fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony") 
Fibro <- FindNeighbors(Fibro, reduction = "harmony") %>% FindClusters() %>% RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony") 
DimPlot(Fibro, group.by = "project", reduction = "umap.harmony") + DimPlot(Fibro, group.by = "seurat_clusters", reduction = "umap.harmony", label = T) 
FacetDotPlot(Fibro, markers_for_stromal, group_by = "seurat_clusters") 
FacetDotPlot(Fibro, kidney_markers, group_by = "seurat_clusters") 

Fibro_filter <- Fibro[, Fibro$seurat_clusters %in% c("0", "1", "4", "6","8")] # Doublet-filtered
FacetDotPlot(Fibro_filter, kidney_markers, group_by = "seurat_clusters") 
Fibro_filter <- NormalizeData(Fibro_filter) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
ElbowPlot(Fibro, ndims = 50) 
Fibro_filter <- IntegrateLayers(Fibro_filter, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony") 
Fibro_filter <- FindNeighbors(Fibro_filter, reduction = "harmony") %>% FindClusters() %>% RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony") 
DimPlot(Fibro_filter, group.by = "project", reduction = "umap.harmony") 

DimPlot(Fibro_filter, reduction = "umap.harmony", group.by = "seurat_clusters", label = T) 
new.cluster.id <- setNames(rep(NA, 9), 0:8) 
new.cluster.id[c("1", "2", "3", "4")] <- "Peri" 
new.cluster.id[c("5")] <- "Fibro" 
new.cluster.id[c("4", "6")] <- "VDR+ Peri" 
new.cluster.id[c("8")] <- "Myo" 
new.cluster.id[c("0", "7")] <- "Peri-Myo transition" 
Fibro_filter$cell_subtype2 <- unname(new.cluster.id[as.character(Fibro_filter$seurat_clusters)]) 
Fibro_filter$cell_subtype2 <- factor(Fibro_filter$cell_subtype2, levels = c("Peri", "VDR+ Peri", "Peri-Myo transition", "Myo", "Fibro")) 

p1 <- FacetDotPlot(Fibro_filter, markers_for_stromal, group_by = "cell_subtype2") 
gs(p1, fig_fig4,  "Fig4G_public_stromal_marker_dotplot", width = 7.5, height = 3.8) 

colors <- c("Peri" = "#9dd279", "Peri-Myo transition" = "#fdcd8b", "VDR+ Peri" = "#ffed9e", "Myo" = "#e19d99", "Fibro" = "#99badd") 
p <- plot_scanpy_umap(Fibro_filter, group_by = "cell_subtype2", label = F, reduction = "umap.harmony",custom_colors = colors) 
gs(p, fig_fig4, "Fig4F_public_stromal_UMAP_legend", width = 6, height = 6) legend_only <- get_legend(p) 
gs(legend_only, fig_fig4, "Fig4F_public_stromal_UMAP_legend", width = 5, height = 6) 

p <- plot_scanpy_umap(Fibro_filter, group_by = "project", label = F, reduction = "umap.harmony", custom_colors = c("#f8766d", "#7cae00", "#00bfc4", "#c77cff")) 
gs(p, fig_supp, "SFig6D_public_stromal_batch_UMAP", 6, 6)
legend_only <- get_legend(p) 
gs(legend_only, fig_supp,  "SFig6D_public_stromal_batch_UMAP_legend", 5, 6)

saveRDS(Fibro_filter, file.path(data_dir, "public_Fibro_filter.rds"))
## ---------- ECM score ---------- 
ECMSigs <- toupper(c('Thbs2', 'Col4a4', 'Tnr', 'Col6a6', 'Vwf', 'Gp6', 'Col1a2', 'Itga3', 'Itga9', 'Col4a5', 'Gp5', 'Lamc1', 'Itga6', 'Itga1', 'Itgb5', 'Itgb8', 'Sv2c', 'Thbs3', 'Hmmr', 'Cd44', 'Lama3', 'Lamb1', 'Lamb3', 'Col6a4', 'Lama4', 'Itgb4', 'Lamc3', 'Itgb7', 'Col1a1', 'Fn1', 'Vtn', 'Sv2a', 'Col9a2', 'Itgb3', 'Gp1bb', 'Col9a3', 'Hspg2', 'Npnt', 'Gp1ba', 'Col9a1', 'Itga5', 'Itgav', 'Lama2', 'Col2a1', 'Spp1', 'Sdc4', 'Tnxb', 'Itga10', 'Col4a3', 'Sdc1', 'Lama5', 'Sv2b', 'Col6a3', 'Reln', 'Itga7', 'Itga2', 'Lamc2', 'Col6a2', 'Ibsp', 'Col4a6', 'Comp', 'Col6a1', 'Tnn', 'Itga11', 'Itgb6', 'Lamb2', 'Chad', 'Thbs4', 'Gp9', 'Itga8', 'Dag1', 'Itga2b', 'Agrn', 'Col4a2', 'Col4a1', 'Cd36', 'Lama1', 'Col6a5', 'Itgb1', 'Itga4', 'Thbs1', 'Tnc', 'Cd47'))
ECMSigs <- ECMSigs[ECMSigs %in% rownames(Fibro_filter)]
ECMSigs <- list('ECM' = ECMSigs)
Fibro_filter <- AddModuleScore_UCell(Fibro_filter, assay='RNA',slot = 'data',features = ECMSigs, ncores = 40, maxRank=2000)

## ---------- CytoTRACE score ---------- 
Fibro_filter <- JoinLayers(Fibro_filter)
mat <- as.matrix(GetAssayData(Fibro_filter, assay = "RNA", layer = "counts"))
CytoTRACE_score <- CytoTRACE(mat = mat, enableFast = F, ncores = 20)
Fibro_filter$CytoTRACE_score <- CytoTRACE_score$CytoTRACE[colnames(Fibro_filter)]

df <- Fibro_filter@meta.data 
df$cell_subtype2 <- factor(df$cell_subtype2,
                           levels = c("Myo", "Fibro", "Peri-Myo transition", "Peri", "VDR+ Peri"))
celltypes <- as.character(unique(df$cell_subtype2))
comparisons <- combn(celltypes, 2, simplify = FALSE)
comparisons <- comparisons[c(10,3,6,8)]

p <- df %>%
  ggplot(aes(x = cell_subtype2, y = CytoTRACE_score ,fill = cell_subtype2)) +
  geom_violin(trim = F, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_fill_manual(values = colors)+
  theme(text = element_text(size = 14),
        axis.text = element_text(color = "black"),
        plot.caption = element_text(size = 8),
        legend.position = "none",
        axis.line = element_line(colour = "black",linewidth = rel(1)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank()) +
  geom_signif(comparisons = comparisons,
              test = "t.test",
              step_increase=0.1, # test.args = c("greater"), # 单侧检验 
              map_signif_level = c(
                "*" = 0.05,
                "**" = 0.01,
                "***" = 0.001,
                "****" = 0.0001
              )
  )+
  ylab('CytoTRACE score') +
  xlab("")

gs(p, pic_fig4, "Fig4H_public_stromal_CytoTRACE_score", width = 7, height = 7)
saveRDS(Fibro_filter, file.path(data_dir, "public_Fibro_filter_umap_ECM.rds"))

## ---------- Monocle2 ---------- 
monocds <- RunMonocle(Fibro_filter)
monocds <- reduceDimension(
  monocds, max_components = 2, residualModelFormulaStr = "~project", 
  auto_param_selection = FALSE, ncenter = 150, reduction_method = "DDRTree")
monocds <- orderCells(monocds)
monocds <- orderCells(monocds, reverse = T)
saveRDS(monocds, file.path(data_dir, "public_monocds.rds"))

p1 <- plot_cell_trajectory(monocds, color_by = "Pseudotime", cell_size = 0.25) + 
  scale_color_gradientn(colors = colorRampPalette(c("#1e3a7e", "white", "#d81d1e"))(200)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
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
  labs(color = "Pseudotime", title = "")
p1
gs(p1, pic_fig4, "Fig4J_pseudotime", width = 5.45, height = 4.5)

p2 <- plot_cell_trajectory(monocds, color_by = "cell_subtype2", cell_size = 0.45) + 
  scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
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
  labs(color = "Cell Subtype", title = "")
gs(p2, pic_fig3, "Fig4I_monocle2_cellsubtype", width = 6, height = 4.5)

df <- data.frame(
  Pseudotime = monocds$Pseudotime,
  ECM_Ucell  = monocds$ECM_UCell,
  celltype = monocds$cell_subtype2,
  State = monocds$State
)
plot_cell_trajectory(monocds, color_by = "State", cell_size = 0.45)
df <- df %>% filter(State %in% c(1,7))
p2 <- ggplot(df, aes(Pseudotime, ECM_Ucell, color = celltype)) +
  geom_point(alpha = 1, size = 0.65) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
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
  labs(color = "Cell Subtype", title = "")
p2
gs(p2, pic_fig4, "Fig4K_ECMscore", width = 5, height = 3.8)

# trajectory genes
Time_diff <- differentialGeneTest(monocds, cores = 20,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- filter(Time_diff,qval < 0.05)
saveRDS(Time_diff, "03_analysis/Trajectory/public_stromal_time_diff.rds")

Time_diff <- readRDS("03_analysis/Trajectory/public_stromal_time_diff.rds")
Time_genes <- top_n(Time_diff, n = 100, -qval) %>% pull(gene_short_name) %>% as.character()

ecm_genes <- c(
  "COL1A1","COL1A2","COL3A1","COL6A2","COL6A3","COL14A1", "FN1","LUM","DCN","VCAN","FBLN1","FBLN5","MFAP4","MGP",
  "TIMP1","PCOLCE","MMP2","CTGF","IGFBP7","IGFBP4","IGFBP6","IGF2","SFRP4","VIM","SRPX")
p <- plot_pseudotime_heatmap(monocds[ecm_genes,],return_heatmap = T,
                        num_clusters = 1,
                        cluster_rows = T,
                        hmcols = colorRampPalette(c("#1e3a7e", "white", "#d81d1e"))(200),
                        cores = 1,
                        show_rownames = T) 
gs(p, fig_supp, "SFig6E_trajectory_genes(human)", width = 5, height = 7)
