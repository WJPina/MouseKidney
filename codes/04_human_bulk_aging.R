#!/usr/bin/env Rscript
source("common.R")
suppressPackageStartupMessages(library(readxl))
df <- as.data.frame(read_excel(input("bulk", "VDR_CDKN1A_AGE.xlsx")))[, c(1, 2, 4)]
names(df) <- c("age", "CDKN1A", "VDR")
df[] <- lapply(df, as.numeric)
df <- df[complete.cases(df), ]
traj <- do.call(rbind, lapply(c("VDR", "CDKN1A"), function(g) {
  med <- aggregate(df[[g]], list(age = df$age), median); names(med)[2] <- "value"
  med$z <- as.numeric(scale(med$value))
  age <- seq(min(med$age), max(med$age), by = 1)
  data.frame(age = age, gene = g, z = as.numeric(scale(predict(loess(z ~ age, med, span = 0.7), data.frame(age = age)))))
}))
save_table(df, "Supp_human_aging_samples"); save_table(traj, "Supp_human_aging_trajectories")
save_plot(ggplot(traj, aes(age, z, color = gene)) + geom_line() + theme_classic(), "Supp_human_aging_trajectories")
ct <- cor.test(df$CDKN1A, df$VDR, method = "spearman", exact = FALSE)
save_table(data.frame(rho = unname(ct$estimate), p_value = ct$p.value, n = nrow(df)), "Supp_human_aging_correlation")
save_plot(ggplot(df, aes(CDKN1A, VDR)) + geom_point() + geom_smooth(method = "lm") + theme_classic() +
  labs(subtitle = sprintf("Spearman rho = %.3f, P = %.3g", ct$estimate, ct$p.value)), "Supp_human_aging_correlation")
