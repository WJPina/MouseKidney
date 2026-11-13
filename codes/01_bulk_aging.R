#!/usr/bin/env Rscript
source("common.R")
suppressPackageStartupMessages(library(DESeq2))
counts <- read.csv(input("bulk", "GSE132040_counts.csv.gz"), row.names = 1, check.names = FALSE)
colnames(counts) <- sub("\\.gencode\\.vM19", "", colnames(counts))
counts <- as.matrix(counts[!grepl("__", rownames(counts)), ])
md <- read.csv(input("bulk", "GSE132040_metadata.csv.gz"), check.names = FALSE)
rownames(md) <- md[["Sample name"]]
md$age <- suppressWarnings(as.numeric(md[["characteristics: age"]]))
md$tissue <- gsub("_[0-9]+", "", md[["source name"]])
md <- md[is.finite(md$age) & !grepl("^NA", md$tissue), ]
samples <- intersect(rownames(md), colnames(counts))
counts <- counts[rowSums(counts[, samples] > 0) >= 10, samples]
counts <- counts[, colSums(counts) > quantile(colSums(counts), 0.05)]
md <- md[colnames(counts), ]
dds <- DESeqDataSetFromMatrix(counts, md, design = ~1)
mat <- assay(vst(dds, blind = TRUE))
md$Vdr <- mat["Vdr", ]; md$Cdkn1a <- mat["Cdkn1a", ]
save_table(md, "Fig1_normalized_samples", TRUE)
save_plot(ggplot(subset(md, tissue != "WBC"), aes(reorder(tissue, -Vdr, median), Vdr)) +
  geom_boxplot() + theme_classic() + labs(x = NULL, y = "VST expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)), "Fig1A_Vdr_tissues", 9, 4)

# Preserve the archived median-by-age, span=0.7 LOESS calculation.
trajectory <- function(df, gene, by = 0.25) {
  med <- aggregate(df[[gene]], list(age = df$age), median)
  names(med)[2] <- "expression"
  med$z <- as.numeric(scale(med$expression))
  age <- seq(min(med$age), max(med$age), by = by)
  z <- predict(loess(z ~ age, med, span = 0.7), data.frame(age = age))
  data.frame(age = age, gene = gene, z = as.numeric(scale(z)))
}
traj <- do.call(rbind, lapply(c("Kidney", "Skin", "Small_Intestine"), function(tissue) {
  z <- trajectory(md[md$tissue == tissue, ], "Vdr")
  z$tissue <- tissue; z
}))
save_table(traj, "Fig1B_Vdr_age_trajectories")
save_plot(ggplot(traj, aes(age, z, color = tissue)) + geom_line() + theme_classic(), "Fig1B_Vdr_age_trajectories")
kidney <- md[md$tissue == "Kidney", ]
traj <- rbind(trajectory(kidney, "Vdr"), trajectory(kidney, "Cdkn1a"))
save_table(traj, "Fig1C_kidney_trajectories")
save_plot(ggplot(traj, aes(age, z, color = gene)) + geom_line() + theme_classic(), "Fig1C_kidney_trajectories")
ct <- cor.test(kidney$Cdkn1a, kidney$Vdr, method = "spearman", exact = FALSE)
save_table(data.frame(rho = unname(ct$estimate), p_value = ct$p.value, n = nrow(kidney)), "Fig1D_correlation")
save_plot(ggplot(kidney, aes(Cdkn1a, Vdr)) + geom_point() + geom_smooth(method = "lm") +
  theme_classic() + labs(subtitle = sprintf("Spearman rho = %.3f, P = %.3g", ct$estimate, ct$p.value)), "Fig1D_correlation")
