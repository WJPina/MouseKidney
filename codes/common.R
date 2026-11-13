suppressPackageStartupMessages(library(ggplot2))
data_dir <- "../data"
out_dir <- "../results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(12135)
cores <- 4

input <- function(...) file.path(data_dir, ...)
save_table <- function(x, name, row.names = FALSE) {
  if (row.names) x <- cbind(feature = rownames(x), x)
  write.csv(x, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
}
save_plot <- function(p, name, width = 6, height = 4) {
  for (ext in c("pdf", "png")) ggsave(file.path(out_dir, paste0(name, ".", ext)),
    p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
}
rna <- function(x, layer = "data") SeuratObject::GetAssayData(x, assay = "RNA", layer = layer)
read_object <- function(...) {
  x <- readRDS(input(...))
  DefaultAssay(x) <- "RNA"
  JoinLayers(x)
}
condition_colors <- c(WT = "#9fa0be", KO = "#e68f8b")
stromal_colors <- c(Peri = "#fbdd9a", Myo = "#e19d99", Fibro = "#99badd",
  "VDR+ Peri" = "#ffed9e", "Peri-Myo transition" = "#fdcd8b")
stromal_markers <- list(Peri = c("Pdgfrb", "Rgs5", "Notch3", "Mylk"),
  Myo = c("Acta2", "Tagln", "Col1a2", "Fn1"), Fibro = c("Col1a1", "Col3a1", "Dcn", "Lum"))
violin <- function(md, group, value, colors = NULL) {
  p <- ggplot(md, aes(x = .data[[group]], y = .data[[value]], fill = .data[[group]])) +
    geom_violin(trim = FALSE, linewidth = 0) + geom_boxplot(width = 0.12, outlier.shape = NA) +
    theme_classic() + labs(x = NULL, y = value) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  if (!is.null(colors)) p <- p + scale_fill_manual(values = colors)
  p
}
ttest <- function(a, b, alternative = "two.sided") {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (min(length(a), length(b)) < 2 || (sd(a) == 0 && sd(b) == 0)) return(NA_real_)
  t.test(a, b, alternative = alternative)$p.value
}
go_enrichment <- function(genes, background, name) {
  ids <- function(g) unique(na.omit(AnnotationDbi::mapIds(org.Mm.eg.db, keys = unique(g),
    keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")))
  if (!length(genes)) return(invisible(NULL))
  tab <- as.data.frame(enrichGO(ids(genes), OrgDb = org.Mm.eg.db, universe = ids(background), ont = "ALL"))
  save_table(tab, name)
  tab <- tab[!is.na(tab$qvalue) & tab$qvalue < 0.05, , drop = FALSE]
  tab <- head(tab[order(tab$qvalue), ], 10)
  if (nrow(tab)) save_plot(ggplot(tab, aes(-log10(qvalue), reorder(Description, -qvalue), size = Count)) +
    geom_point(color = "#d81d1e") + theme_classic() + labs(y = NULL), name, 9, 5)
}
raw_stromal <- function(x) {
  atlas <- read_object("mouse", "kidney.rds")
  counts <- rna(atlas, "counts")[intersect(rownames(x), rownames(atlas)), colnames(x)]
  y <- NormalizeData(CreateSeuratObject(counts, meta.data = x[[]]), verbose = FALSE)
  for (red in names(x@reductions)) y[[red]] <- x[[red]]
  y
}
build_trajectory <- function(x, group, origin) {
  counts <- rna(x, "counts")
  pd <- Biobase::AnnotatedDataFrame(x[[]])
  fd <- Biobase::AnnotatedDataFrame(data.frame(gene_short_name = rownames(counts), row.names = rownames(counts)))
  cds <- newCellDataSet(as(counts, "sparseMatrix"), phenoData = pd, featureData = fd,
    lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
  cds <- detectGenes(estimateDispersions(estimateSizeFactors(cds)), min_expr = 0.1)
  cds <- cds[Biobase::fData(cds)$num_cells_expressed >= 10, ]
  disp <- dispersionTable(cds)
  ordering <- with(disp, gene_id[mean_expression >= 0.05 & dispersion_empirical >= dispersion_fit])
  cds <- setOrderingFilter(cds, ordering)
  cds <- orderCells(reduceDimension(cds, max_components = 2, reduction_method = "DDRTree"))
  tab <- table(Biobase::pData(cds)$State, Biobase::pData(cds)[[group]] == origin)
  state <- as.numeric(rownames(tab)[which.max(tab[, "TRUE"] / rowSums(tab))])
  orderCells(cds, root_state = state)
}
plot_trajectory <- function(cds, group, prefix, branch_states = NULL) {
  save_plot(plot_cell_trajectory(cds, color_by = group, cell_size = 0.5) +
    scale_color_manual(values = stromal_colors, name = "Cell subtype") +
    guides(color = guide_legend(nrow = if (length(unique(Biobase::pData(cds)[[group]])) > 3) 2 else 1)),
    paste0(prefix, "_celltypes"))
  save_plot(plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5), paste0(prefix, "_pseudotime"))
  md <- Biobase::pData(cds)
  save_table(md, paste0(prefix, "_cells"), TRUE)
  if (!is.null(branch_states)) md <- md[md$State %in% branch_states, ]
  save_plot(ggplot(md, aes(Pseudotime, ECM_UCell, color = .data[[group]])) +
    geom_point(size = 0.6) + geom_smooth(method = "loess", color = "red", se = FALSE) +
    scale_color_manual(values = stromal_colors) + theme_classic(), paste0(prefix, "_ECM"))
}
trajectory_genes <- function(cds, genes, prefix) {
  de <- differentialGeneTest(cds, cores = cores, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  save_table(de, paste0(prefix, "_gene_tests"), TRUE)
  pdf(file.path(out_dir, paste0(prefix, "_heatmap.pdf")), width = 5, height = 7)
  print(plot_pseudotime_heatmap(cds[intersect(genes, rownames(cds)), ], return_heatmap = TRUE,
    num_clusters = 1, cluster_rows = TRUE, cores = 1, show_rownames = TRUE))
  dev.off()
}
