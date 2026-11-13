#!/usr/bin/env Rscript
source("common.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(mclust)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# Atlas, expression scores and cell composition.
sce <- read_object("mouse", "kidney.rds")
sce$Condition <- factor(sce$Condition, c("WT", "KO"))
Idents(sce) <- "celltype"
save_plot(DimPlot(sce, group.by = "celltype", raster = TRUE), "Fig3B_mouse_atlas", 8, 6)
save_table(as.data.frame(Embeddings(sce, "umap")), "Fig3B_umap", TRUE)
markers <- list(PT = c("Slc27a2", "Lrp2"), DCT = c("Pgam2", "Pvalb", "Wnk1"),
  DLH = c("Epha7", "Bst1"), ALH = c("Tmem207", "Slc12a1", "Umod"),
  CD_IC = c("Uqcrb", "Slc25a4", "Car2", "Atp6v1g3", "Atp6v0d2"),
  CD_PC = c("Aqp3", "Apela", "Fxyd4", "Hsd11b2", "Aqp2"), Podo = c("Cdkn1c", "Magi2", "Nphs2", "Mafb"),
  Endo = c("Tek", "Adgrl4", "Emcn", "Plvap"), Stromal = c("Rgs5", "Pdgfrb", "Acta2", "Col1a1"),
  UE = c("Ly6d", "Upk3a", "Psca"), RBC = c("Hba-a2", "Hba-a1", "Hbb-bs"))
p <- DotPlot(sce, features = markers)
save_plot(p + RotatedAxis(), "Fig3C_markers", 13, 4)
save_table(p$data, "Fig3C_markers")
genes <- c("Cdkn1a", "Trp53", "Tgfb1", "Fn1", "Acta2", "Vim")
expr <- as.data.frame(t(as.matrix(rna(sce)[genes, ])))
expr$Condition <- sce$Condition
save_table(data.frame(gene = genes, p_value = vapply(genes, function(g)
  ttest(expr[expr$Condition == "KO", g], expr[expr$Condition == "WT", g], "greater"), numeric(1))), "Fig3D_cell_tests")
save_plot(wrap_plots(lapply(genes, function(g) violin(expr, "Condition", g, condition_colors)), nrow = 1), "Fig3D_senescence_genes", 14, 4)
save_plot(violin(sce[[]], "Condition", "ECM_UCell", condition_colors), "Fig3F_ECM")
save_plot(FeaturePlot(sce, "ECM_UCell", order = TRUE, min.cutoff = 0.05, max.cutoff = 0.225), "Fig3G_ECM_umap", 6, 5)
save_plot(violin(sce[[]], "Condition", "Senescence_UCell", condition_colors), "Supp_senescence_score")

# Rank components by their fitted means before applying the published 1-4/5-7/8 grouping.
gmm <- Mclust(sce$ECM_UCell, G = 8, verbose = FALSE)
component <- match(gmm$classification, order(as.numeric(gmm$parameters$mean)))
sce$ECM_group <- factor(ifelse(component <= 4, "low", ifelse(component <= 7, "middle", "high")), c("low", "middle", "high"))
saveRDS(gmm, file.path(out_dir, "Fig3H_gmm.rds"))
save_plot(ggplot(sce[[]], aes(ECM_UCell, fill = ECM_group)) + geom_histogram(bins = 60) + theme_classic(), "Fig3H_ECM_groups")
tab <- as.data.frame(prop.table(table(sce$celltype_fine, sce$ECM_group), 1))
names(tab) <- c("celltype", "ECM_group", "fraction")
save_table(tab, "Fig3I_ECM_composition")
save_plot(ggplot(subset(tab, !celltype %in% c("RBC", "Unknow")), aes(celltype, fraction, fill = ECM_group)) +
  geom_col() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)), "Fig3I_ECM_composition", 10, 4)

composition <- function(x, group, name) {
  tab <- as.data.frame(prop.table(table(x$orig.ident, x[[group, drop = TRUE]]), 1))
  names(tab) <- c("sample", "celltype", "fraction")
  tab$Condition <- substr(as.character(tab$sample), 1, 2)
  tests <- do.call(rbind, lapply(split(tab, tab$celltype), function(d)
    data.frame(celltype = as.character(d$celltype[1]), p_value = ttest(d$fraction[d$Condition == "KO"], d$fraction[d$Condition == "WT"]))))
  tests$p_adjust_BH <- p.adjust(tests$p_value, "BH")
  save_table(tab, name); save_table(tests, paste0(name, "_tests"))
  save_plot(ggplot(tab, aes(celltype, fraction, color = Condition)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.5)) +
    stat_summary(aes(group = Condition), fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
    scale_color_manual(values = condition_colors) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)), name, 9, 4)
}
composition(sce, "celltype", "Fig3J_broad_proportions")
composition(read_object("mouse", "stromal.rds"), "celltype", "Fig3J_stromal_proportions")
save_plot(VlnPlot(sce, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", pt.size = 0), "Supp_QC_after", 11, 4)
save_plot(DimPlot(sce, group.by = "orig.ident", raster = TRUE), "Supp_mouse_batch", 7, 5)
save_plot(DimPlot(sce, group.by = "celltype_fine", raster = TRUE), "Supp_mouse_fine_atlas", 9, 6)
save_table(sce[[]], "mouse_metadata", TRUE)

# Differential expression and GO enrichment.
# The supplement specifies cell-level Wilcoxon tests with Bonferroni correction.
Idents(sce) <- "Condition"
de <- FindMarkers(sce, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0.1, min.pct = 0.01)
save_table(de, "Fig3_all_Wilcoxon_DE", TRUE)
go_enrichment(rownames(de)[de$p_val_adj < 0.05 & de$avg_log2FC > 1], rownames(sce), "Fig3E_GO_up")
go_enrichment(rownames(de)[de$p_val_adj < 0.05 & de$avg_log2FC < -1], rownames(sce), "Supp_GO_down")

st <- read_object("mouse", "stromal.rds")
st$group <- paste(st$celltype, st$Condition)
Idents(st) <- "group"
for (ct in c("Peri", "Myo", "Fibro")) {
  de <- FindMarkers(st, ident.1 = paste(ct, "KO"), ident.2 = paste(ct, "WT"), min.pct = 0.25)
  save_table(de, paste0("Fig3K_", ct, "_DE"), TRUE)
  de$gene <- rownames(de)
  de$direction <- ifelse(de$p_val < 0.05 & abs(de$avg_log2FC) > 1, ifelse(de$avg_log2FC > 0, "KO", "WT"), "NS")
  save_plot(ggplot(de, aes(avg_log2FC, -log10(pmax(p_val, 1e-300)), color = direction)) +
    geom_point(size = 0.6) + theme_classic() + labs(title = ct), paste0("Fig3K_", ct, "_volcano"))
  go_enrichment(de$gene[de$p_val < 0.05 & de$avg_log2FC > 1], rownames(st), paste0("Supp_", ct, "_GO_up"))
}

# Proximal-tubule and endothelial subtypes.
markers <- list(PT = list(S1 = c("Slc5a2", "Slc5a12", "Spp2"),
  S2 = c("Slc22a8", "Slco1a6", "Slc13a3", "Slc22a6"), S3 = c("Slc7a13", "Slc27a2", "Slc7a12", "Slc22a7")),
  Endo = list(Endo = c("Pecam1", "Cdh5", "Kdr"), Arteriole = c("Hey2", "Gata6", "Gper1", "Adh7"),
    Venule = c("Aqp1", "Lef1", "Slc14a1", "Hpgd"), Capillary = c("Plvap", "Esm1", "Col15a1", "Vwa1")))
for (ct in names(markers)) {
  x <- sce[, sce$celltype == ct & !sce$celltype_fine %in% c("Unknow", "Unknown")]
  x[["RNA"]] <- split(x[["RNA"]], f = x$orig.ident)
  x <- RunPCA(ScaleData(FindVariableFeatures(NormalizeData(x))), verbose = FALSE)
  x <- IntegrateLayers(x, method = CCAIntegration, new.reduction = "cca", verbose = FALSE)
  x <- FindNeighbors(x, reduction = "cca", dims = 1:20, verbose = FALSE)
  x <- FindClusters(x, resolution = 0.8, random.seed = 12135, verbose = FALSE)
  x <- RunUMAP(x, reduction = "cca", dims = 1:20, seed.use = 12135, verbose = FALSE)
  x <- JoinLayers(x)
  # Use the released cell-level annotations for the refitted embeddings.
  Idents(x) <- "celltype_fine"
  save_plot(DimPlot(x, group.by = "celltype_fine", raster = TRUE), paste0("Supp_", ct, "_subtypes"))
  save_plot(DotPlot(x, features = markers[[ct]]) + RotatedAxis(), paste0("Supp_", ct, "_markers"), 8, 4)
  save_plot(DimPlot(x, group.by = "orig.ident", raster = TRUE), paste0("Supp_", ct, "_batch"))
}
