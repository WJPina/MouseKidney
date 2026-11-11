## ===============================
## Figure: Fig3 – PT cell subtypes
## ===============================

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

## ---------- Paths ----------
fig_dir   <- here("figures")
fig_fig3  <- here("figures", "Fig3")
fig_supp  <- here("figures", "Supplementary")
data_dir  <- here("data", "processed", "scRNA")

dir.create(fig_fig3, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp, recursive = TRUE, showWarnings = FALSE)

## ---------- Source custom functions ----------
source(here("R", "sc_yuanplot.R"))

## ---------- Save function ----------
gs <- function(plot, path, name, width = 5, height = 5){
  ggsave(file.path(path, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 600)
  ggsave(file.path(path, paste0(name, ".pdf")),
         plot, width = width, height = height)
}

## ---------- Load data ----------
sce <- readRDS(file.path(data_dir, "inhouse_sce_filter.rds"))

## ---------- PT extraction & preprocessing ----------
PT <- sce[, sce$celltype == "PT"]

PT <- PT |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA()

ElbowPlot(PT, ndims = 50)

PT <- PT |>
  FindNeighbors(dims = 1:20) |>
  FindClusters(dims = 1:20) |>
  RunUMAP(dims = 1:20)

## ---------- Marker visualization ----------
KidneyMarkers <- list(
  PT   = c("Slc27a2", "Lrp2"),
  DCT  = c("Pgam2", "Pvalb", "Wnk1"),
  DLH  = c("Epha7", "Bst1"),
  ALH  = c("Tmem207", "Slc12a1", "Umod"),
  CD_IC = c("Uqcrb", "Slc25a4", "Car2", "Atp6v1g3", "Atp6v0d2"),
  CD_PC = c("Aqp3", "Apela", "Fxyd4", "Hsd11b2", "Aqp2"),
  Podo = c("Cdkn1c", "Magi2", "Nphs2", "Mafb"),
  Endo = c("Tek", "Adgrl4", "Emcn", "Plvap"),
  Stromal = c("Rgs5", "Pdgfrb", "Acta2", "Col1a1"),
  UE   = c("Ly6d", "Upk3a", "Psca"),
  RBC  = c("Hba-a2", "Hba-a1", "Hbb-bs")
)

FacetDotPlot(PT, KidneyMarkers, group_by = "seurat_clusters")
DimPlot(PT, group.by = "orig.ident")

## ---------- Batch correction ----------
PT[["RNA"]] <- split(PT[["RNA"]], f = PT$orig.ident)

PT <- PT |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA()

PT <- IntegrateLayers(
  PT,
  method = CCAIntegration,
  new.reduction = "CCA"
)

PT <- PT |>
  FindNeighbors(reduction = "CCA", dims = 1:20) |>
  FindClusters(dims = 1:20) |>
  RunUMAP(
    dims = 1:20,
    reduction = "CCA",
    reduction.name = "umap.cca"
  )

## ---------- PT subtype annotation ----------
pt_markers <- list(
  S1  = c("Slc5a2", "Slc5a12", "Spp2"),
  S2  = c("Slc22a8", "Slco1a6", "Slc13a3", "Slc22a6"),
  S3  = c("Slc7a13", "Slc27a2", "Slc7a12", "Slc22a7")
)

FacetDotPlot(PT, pt_markers, group_by = "seurat_clusters")

new.cluster.id <- setNames(
  rep("PT-S2", length(unique(PT$seurat_clusters))),
  unique(PT$seurat_clusters)
)

new.cluster.id[c("8", "10", "11", "14", "15")] <- "PT-S1"
new.cluster.id["13"] <- "PT-S3"

PT$cellsubtype <- unname(new.cluster.id[as.character(PT$seurat_clusters)])

## ---------- Plot Fig3 ----------
pt_colors <- c(
  "PT-S1" = "#d0ead5",
  "PT-S2" = "#e6c4df",
  "PT-S3" = "#c2c0e6"
)

p_fig3 <- plot_scanpy_umap(
  PT,
  group_by = "cellsubtype",
  reduction = "umap.cca",
  custom_colors = pt_colors
)

gs(p_fig3, fig_fig3, "Fig3B_PT_subtype_UMAP", width = 5.5, height = 6)

## ---------- Supplementary ----------
p_dot <- FacetDotPlot(
  PT,
  pt_markers,
  group_by = "cellsubtype",
  ylab = "Cell subtype"
)

gs(p_dot, fig_supp, "SFig4_PT_marker_dotplot", width = 5.5, height = 3.2)

## ---------- Save results ----------
sce$celltype_fine <- as.character(sce$celltype_fine)
sce$celltype_fine[sce$celltype == "PT"] <- "Unknown"
sce$celltype_fine[colnames(PT)] <- PT$cellsubtype

saveRDS(sce, file.path(data_dir, "inhouse_sce_filter.rds"))
saveRDS(PT,  file.path(data_dir, "inhouse_PT.rds"))
