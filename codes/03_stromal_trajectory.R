#!/usr/bin/env Rscript
source("common.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
})
`%||%` <- function(x, y) if (is.null(x)) y else x

# Mouse stromal cells: stored scores and a trajectory refitted from RNA UMIs.
x <- read_object("mouse", "stromal.rds")
Idents(x) <- "celltype"
save_plot(DimPlot(x, group.by = "celltype", cols = stromal_colors[c("Peri", "Myo", "Fibro")]), "Supp_stromal_umap")
save_plot(DotPlot(x, features = stromal_markers) + RotatedAxis(), "Supp_stromal_markers", 8, 4)
save_plot(violin(x[[]], "celltype", "CytoTRACE_score", stromal_colors), "Fig4A_CytoTRACE")
save_table(x[[]], "Fig4A_CytoTRACE", TRUE)
pairs <- list(c("Peri", "Myo"), c("Peri", "Fibro"), c("Myo", "Fibro"))
save_table(do.call(rbind, lapply(pairs, function(z) data.frame(group1 = z[1], group2 = z[2],
  p_value = ttest(x$CytoTRACE_score[x$celltype == z[1]], x$CytoTRACE_score[x$celltype == z[2]], "greater")))), "Fig4A_tests")
cds <- build_trajectory(raw_stromal(x), "celltype", "Peri")
saveRDS(cds, file.path(out_dir, "mouse_monocle.rds"))
plot_trajectory(cds, "celltype", "Fig4BCD_mouse")
trajectory_genes(cds, c("Col3a1", "Col5a2", "Fn1", "Fbn1", "Timp1", "Bgn", "Dcn", "Olfml3", "Col15a1",
  "Tnc", "Spp1", "Serpine1", "Ctgf", "Col1a2", "Fbln2", "Mfap5", "Emp1", "Mgp", "Tgfb1", "Cx3cl1",
  "Cdh11", "Loxl1", "Aebp1", "Sparc", "Col1a1", "Cxcl14", "F2r", "Ltbp2", "Fstl1", "Col8a1", "Thbs1", "S100a6"), "Fig4E_mouse")

# Human stromal cells: published cell labels, scores and trajectory coordinates.
x <- read_object("human", "public_Fibro_filter.rds")
cds <- readRDS(input("human", "public_monocds.rds"))
Idents(x) <- "cell_subtype2"
save_plot(DimPlot(x, reduction = "umap.harmony", group.by = "cell_subtype2", cols = stromal_colors), "Fig4F_human_stromal", 7, 5)
save_table(cbind(as.data.frame(Embeddings(x, "umap.harmony")), x[[]][, c("cell_subtype2", "project")]), "Fig4F_umap", TRUE)
markers <- c(list(VDR = "VDR"), lapply(stromal_markers, toupper))
markers$Myo <- c("ACTA2", "TAGLN", "COL1A2", "POSTN")
p <- DotPlot(x, features = markers) + RotatedAxis()
save_plot(p, "Fig4G_human_markers", 9, 4)
save_table(p$data, "Fig4G_human_markers")
save_plot(violin(x[[]], "cell_subtype2", "CytoTRACE_score", stromal_colors), "Fig4H_human_CytoTRACE", 8, 5)
md <- x[[]]
save_table(md, "human_stromal_metadata", TRUE)
save_table(do.call(rbind, lapply(setdiff(unique(as.character(md$cell_subtype2)), "VDR+ Peri"), function(ct)
  data.frame(group1 = "VDR+ Peri", group2 = ct,
    p_value = ttest(md$CytoTRACE_score[md$cell_subtype2 == "VDR+ Peri"], md$CytoTRACE_score[md$cell_subtype2 == ct], "greater")))), "Fig4H_tests")
plot_trajectory(cds, "cell_subtype2", "Fig4IJK_human", c(1, 7))
save_plot(DimPlot(x, reduction = "umap.harmony", group.by = "project"), "Supp_human_stromal_batch")
