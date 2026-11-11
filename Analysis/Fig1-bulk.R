library(DESeq2)
library(data.table)
library(readr)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
library(patchwork)
library(readxl)
set.seed(12135)

## ---------- Load count matrix ----------
counts <- read_csv("GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
colnames(counts) <- gsub("\\.gencode\\.vM19", "", colnames(counts))
counts <- counts %>% as.data.frame() %>% 
  column_to_rownames(var = "gene") 
counts <- counts[!grepl("__", rownames(counts)), ] 

## ---------- Load metadata ----------
metadata <- read_csv("GSE132040_MACA_Bulk_metadata.csv.gz") %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "Sample name")
metadata$`source name` <- gsub("_[0-9]+", "", metadata$`source name`)
metadata_filter <- metadata[!grepl("NA", metadata$`characteristics: age`) & !grepl("^NA",metadata$`source name`), ]

## ---------- Match samples ----------
common_sample <- intersect(rownames(metadata_filter), colnames(counts))
counts <- counts[, common_sample] 
mdf <- metadata_filter[common_sample, ]

mdf$age <- as.numeric(mdf$`characteristics: age`)
mdf$tissue <- gsub("_[0-9]+", "",mdf$`source name`)

## ---------- Gene and sample filtering ----------
gene_keep <- rowSums(counts > 0) >= 10
counts_filt <- counts[gene_keep, ]
cat("filter", table(gene_keep)[1], "ratained", table(gene_keep)[2])

sample_total <- colSums(counts_filt)
hist(sample_total, breaks = 100)
abline(v = quantile(sample_total, 0.05), col = "red", lwd = 2, lty = 2)
summary(sample_total)
sample_keep <- sample_total > quantile(sample_total, 0.05)
counts_filter <- counts_filt[, sample_keep]
cat("filter", table(sample_keep)[1], "retained", table(sample_keep)[2])

mdf <- mdf[colnames(counts_filter), ]
counts_filter  <- counts_filter[, rownames(mdf)]

## ---------- DESeq2 VST normalization ----------
dds <- DESeqDataSetFromMatrix(
  countData = counts_filter,
  colData = mdf,
  design = ~1
)
vst_mat <- vst(dds, blind = TRUE) %>% assay()
vst_out <- vst_mat %>% as.data.frame() %>% rownames_to_column("Gene")

## ---------- Plot ----------
mdf$Cdkn1a <- as.numeric(vst_mat["Cdkn1a", rownames(mdf)])
mdf$Vdr <- as.numeric(vst_mat["Vdr", rownames(mdf)])

## 1. Fig1A ----
df <- mdf %>% 
  dplyr::select(`source name`, age, Cdkn1a)
df <- df[df$`source name` != "WBC",]
df2 <- df %>% 
  filter(age >= 3)
baseline <- df2 %>%
  filter(age == 3) %>%
  group_by(`source name`) %>%
  summarise(base = mean(Cdkn1a, na.rm = TRUE))

df_plot <- df2 %>%
  left_join(baseline, by = "source name") %>%
  group_by(age, `source name`) %>% 
  summarise(Cdkn1a = mean(Cdkn1a), base = mean(base)) %>% 
  mutate(change = log2(Cdkn1a / base))

df_plot$color_group <- ifelse(df_plot$`source name` == "Kidney", "Kidney", "Other")
df_plot$`source name` <- factor(df_plot$`source name` ,
                                levels = c("BAT", "Bone", "Brain", "GAT", "Heart", "Limb_Muscle", 
                                           "Liver", "Lung", "Marrow", "MAT", "Pancreas", "SCAT", "Skin", 
                                           "Small_Intestine", "Spleen", "Kidney"))
p1 <- ggplot(df_plot, aes(x = age, y = change, color = color_group, group = `source name`)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cdkn1a expression change relative to 3 Months",
    x = "Age (months)",
    y = "log2 fold change",
    color = "Tissue"
  ) +
  scale_x_continuous(breaks = sort(unique(df_plot$age))) +
  scale_color_manual(values = c(
    "Kidney" = "red",
    "Other" = "grey90"
  )) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.1,0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.45)
  )
p1

## 2. Fig1B ----
df <- mdf %>% 
  mutate(source = reorder(`source name`, Vdr, FUN = median, decreasing = T)) %>% 
  dplyr::select(source, Vdr, age)
df <- df[df$source != "WBC",]
p2 <- ggplot(df,aes(x = source, y = Vdr, fill = source)) +
  stat_boxplot(geom = "errorbar", width = 0.45, size = 0.35) +
  geom_boxplot(width = 0.7, size = 0.35, outlier.shape = 21, 
               outlier.size = 1.3, outlier.stroke = 0.35, linewidth = 0.35, fatten = 1.15) + # fatten表示中间均值的线相对于其他线 
  theme_classic(base_size = 14) +
  labs(x = "Tissue", y = "VST-normalized expression", title = "VDR mRNA expression across tissues") +
  scale_fill_manual(values = c("#6ab2c3","#dd6771","#d4efb0", colorRampPalette(brewer.pal(12, "Set3"))(14))) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.45)
  )
p2


## ---------- LOESS ----------
calc_loess_traj <- function(expr, age, span = 0.7,
                            pred_age = seq(min(age), max(age), by = 0.25)) {
  df <- data.frame(age = age, expr = expr)
  med_df <- df %>%
    group_by(age) %>% 
    summarise(med = median(expr, na.rm = TRUE), .groups = "drop") %>%
    arrange(age)
  
  med_df$z <- as.numeric(scale(med_df$med))
  fit <- loess(z ~ age, data = med_df, span = span)
  pred_z <- predict(fit, newdata = data.frame(age = pred_age))
  names(pred_z) <- pred_age
  return(pred_z)
}


tissues_unique <- unique(mdf$tissue)

traj_list <- list()
for(ti in tissues_unique){
  message("Processing tissue: ", ti)
  
  tissue_samples <- rownames(mdf)[mdf$tissue == ti]

  mat <- mdf[tissue_samples, c("Vdr", "Cdkn1a"),drop=FALSE]
  age_vec <- mdf[tissue_samples, "age"]
  
  traj_mat <- sapply(c("Vdr", "Cdkn1a"), function(g){
    calc_loess_traj(mat[, g], age_vec)
  })
  
  rownames(traj_mat) <- pred_ages
  traj_list[[ti]] <- traj_mat
}

Kidney <- traj_list[["Kidney"]][, "Vdr"]
Skin <- traj_list[["Skin"]][, "Vdr"]
Small_Intestine <- traj_list[["Small_Intestine"]][, "Vdr"]
df <- data.frame(
  Kidney = scale(Kidney),
  Skin = scale(Skin),
  Small_Intestine = scale(Small_Intestine)
) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Tissue")

colnames(df)[2:dim(df)[2]] <- paste0("X", colnames(df)[2:dim(df)[2]])
Vdr_trajectories <- df

Vdr_traj_long <- Vdr_trajectories %>% 
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Month",
    values_to = "Zscore",
  ) %>% 
  mutate(Month = as.numeric(gsub("X", "", Month)))

max_kidney <- Vdr_traj_long %>% filter(Tissue == "Kidney") %>% slice_max(Zscore, n =1)
p3 <- ggplot(Vdr_traj_long, aes(x = Month, y = Zscore, color = Tissue)) +
  geom_line(size = 0.45) +
  geom_vline(xintercept = max_kidney$Month, linetype = "dashed",
             color = "#dd6771", linewidth = 0.3)+
  theme_classic(base_size = 12) +
  scale_color_manual(values =  c("Small_Intestine"= "#6ab2c3", "Kidney"  = "#dd6771", "Skin" = "#d4efb0"))+
  theme(
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.ticks = element_line(color = "black",linewidth = 0.25),
    legend.position = c(0.8,0.8),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(
    breaks = c(3, 6, 9, 12, 15, 18, 21, 24),
    labels = c("3", "6", "9", "12", "15", "18", "21","24")
  ) +
  labs(title = "Age-associated Vdr expression in mouse")
p3

Kidney <- traj_list[["Kidney"]][, c("Vdr","Cdkn1a")]
Kidney <- scale(Kidney)
ct <- cor.test(Kidney[,1], Kidney[,2], method = "pearson")
r_value <- round(ct$estimate, 2)
p_text <- ifelse(ct$p.value < 0.05, "P < 0.05",
                 paste0("p = ", signif(ct$p.value, 2)))
annot_text <- paste0("R = ", r_value, "\n", p_text)
vdr_cdkn1a_long <- Kidney %>% as.data.frame() %>% 
  rownames_to_column("Month") %>% 
  mutate(
    Month = as.numeric(Month),
    Tissue = "Kidney",
  ) %>% 
  pivot_longer(
    cols = c(Vdr, Cdkn1a),
    names_to = "Gene",
    values_to = "Zscore"
  )
p4 <- ggplot(vdr_cdkn1a_long, aes(x = Month, y = Zscore, color = Gene)) +
  geom_line()+
  theme_classic(base_size = 12) +
  scale_color_manual(values =  c("Cdkn1a"= "#b2a1cc", "Vdr"  = "#dd6771"))+
  theme(
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.ticks = element_line(color = "black",linewidth = 0.25),
    legend.position = c(0.1,0.8),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(
    breaks = c(3, 6, 9, 12, 15, 18, 21, 24),
    labels = c("3", "6", "9", "12", "15", "18", "21","24")
  ) + 
  annotate("text",
           x = 25, y = max(vdr_cdkn1a_long$Zscore)*0.85,
           label = annot_text, size = 3)+
  labs(title = "Vdr–Cdkn1a expression relationship in mouse")
p4

mouse <- mdf[mdf$`source name` == "Kidney",]
mouse$Age <- mouse$age
mouse$Cdkn1a <- scale(mouse$Cdkn1a)
mouse$Vdr <- scale(mouse$Vdr)

mouse <- mouse[mouse$Vdr != 0,]
mouse_long <- mouse %>%
  pivot_longer(
    cols = c(Cdkn1a, Vdr),
    names_to = "Gene",
    values_to = "Expression"
  )
cor_mouse <- cor.test(
  mouse$Cdkn1a,
  mouse$Vdr,
  method = "spearman",
  use = "complete.obs"
)

r_mouse <- round(cor_mouse$estimate, 2)
p_mouse <- signif(cor_mouse$p.value, 2)
p_text <- ifelse(p_mouse < 0.05, "P < 0.05",
                 paste0("p = ", signif(ct$p.value, 2)))
p6 <- ggplot(mouse, aes(x = Cdkn1a, y = Vdr)) +
  geom_point(size = 2, alpha = 0.7, color = "#4E79A7") +
  geom_smooth(method = "lm", se = TRUE, color = "#4E79A7", linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("R = ", r_mouse, "\n", p_text),
    hjust = 1.1, vjust = 1.2,
    size = 4
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cdkn1a expression (Zscore)",
    y = "Vdr expression (Zscore)",
    title = "Correlation between Cdkn1a and Vdr in mouse kidney"
  ) +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(color = "black",linewidth = 0.6),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
p6

## ---------- human ----------
df <- read_excel("VDR_CDKN1A_AGE.xlsx")
df <- df[,c(1,2,4)] %>% as.data.frame()
colnames(df) <- c("Age", "CDKN1A", "VDR")
Kidney <- sapply(c("VDR", "CDKN1A"), function(g){
  calc_loess_traj(df[, g], as.numeric(df$Age), 
                  pred_age = seq(min(as.numeric(df$Age)),max(as.numeric(df$Age)), by = 1))
})
Kidney <- scale(Kidney)
ct <- cor.test(Kidney[,1], Kidney[,2], method = "pearson")
r_value <- round(ct$estimate, 2)
p_text <- ifelse(ct$p.value < 0.05, "P < 0.05",
                 paste0("p = ", signif(ct$p.value, 2)))
annot_text <- paste0("R = ", r_value, "\n", p_text)
vdr_cdkn1a_long <- Kidney %>% as.data.frame() %>% 
  rownames_to_column("Year") %>% 
  mutate(
    Month = as.numeric(Year),
    Tissue = "Kidney",
  ) %>% 
  pivot_longer(
    cols = c(VDR, CDKN1A),
    names_to = "Gene",
    values_to = "Zscore"
  )
p5 <- ggplot(vdr_cdkn1a_long, aes(x = Month, y = Zscore, color = Gene)) +
  geom_line()+
  theme_classic(base_size = 12) +
  scale_color_manual(values =  c("CDKN1A"= "#b2a1cc", "VDR"  = "#dd6771"))+
  theme(
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    axis.ticks = element_line(color = "black",linewidth = 0.25),
    legend.position = c(0.1,0.8),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(
    breaks = seq(30,90,10),
    labels = seq(30,90,10)
  ) + 
  annotate("text",
           x = 90, y = max(vdr_cdkn1a_long$Zscore)*0.85,
           label = annot_text, size = 3) +
  labs(title = "Vdr–Cdkn1a expression relationship in human", x= "Year")
p5


df <- read_excel("VDR_CDKN1A_AGE.xlsx")
df[2:6] <- scale(df[2:6])
plot_df <- data.frame(CDKN1A = df[[2]], VDR = df[[4]])
cor_human <- cor.test(plot_df$CDKN1A, plot_df$VDR, method = "spearman", use = "complete.obs")

r_human <- round(cor_human$estimate, 2)
p_human <- signif(cor_human$p.value, 2)
p_text <- ifelse(p_human < 0.05, "P < 0.05",
                 paste0("p = ", signif(ct$p.value, 2)))
p7 <- ggplot(plot_df, aes(x = CDKN1A, y = VDR)) +
  geom_point(size = 2, alpha = 0.7, color = "#59A14F") +
  geom_smooth(method = "lm", se = TRUE, color = "#59A14F", linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("R = ", r_human, "\n", p_text),
    hjust = 1.1, vjust = 1.2,
    size = 4
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0(colnames(df2)[2], " expression (Zscore)"),
    y = "VDR expression (Zscore)",
    title = "Correlation between CDKN1A and VDR in human kidney"
  ) +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(color = "black",linewidth = 0.6),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


p <- p1 + p2 + p3 +
  p4 + p6 + p5 +  p7 + 
  plot_layout(design = "
  AAABBBCC
  DDEEFFGG") +
  plot_annotation(tag_levels = "A")
ggsave("Fig1A-G.pdf", p, width = 13, height = 7.5)



