rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pheatmap, tidyverse, RColorBrewer)
load("GSE25066.rda")
DEG <- read.table("GSE25066_diff1.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
stopifnot(identical(rownames(group), colnames(expr)))
table(DEG$group)
sig_genes <- rownames(DEG)[DEG$group != "NOT"]
exp_diff <- expr[sig_genes, ]
pheatmap(exp_diff, annotation_col = group, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = FALSE, border_color = NA)
DEG_filtered <- DEG[DEG$group != "NOT", ]
top_up <- DEG_filtered %>% arrange(desc(log2FC)) %>% head(25)
top_down <- DEG_filtered %>% arrange(log2FC) %>% head(25)
selected_genes <- c(rownames(top_up), rownames(top_down))
exp_diff2 <- expr[selected_genes, ]
annotation_colors <- list(response = c(Responder = "#4DBBD5FF", Non_Responder = "#E64B35FF"))
final_heatmap <- pheatmap(exp_diff2, annotation_col = group, annotation_colors = annotation_colors, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), scale = "row", cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, border_color = NA, fontsize_row = 7, cellwidth = 2, cellheight = 10, gaps_row = 25, main = "Top 50 DEGs", legend = TRUE)
sample_count <- ncol(exp_diff2)
ggsave("DEG_heatmap.pdf", plot = final_heatmap, width = max(6, sample_count * 0.15), height = 20, dpi = 300, units = "cm")
ggsave("DEG_heatmap.tiff", plot = final_heatmap$gtable, width = max(6, sample_count * 0.15), height = 20, dpi = 300, units = "cm", device = "tiff")
