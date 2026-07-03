library(rjags)
library(infercnv)
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)
library(cowplot)
library(Matrix)
output_dir <- "CCND1_Analysis_Output"
dir.create(output_dir, showWarnings = FALSE)
seu_combined <- readRDS("seu_combined_harmony_annotated_renamed-1.rds")
cat("Layers of original object:", Layers(seu_combined, assay = "RNA"), "\n")
seu_combined <- JoinLayers(seu_combined, assay = "RNA")
cat("Layers after merge:", Layers(seu_combined, assay = "RNA"), "\n")
cat("Metadata columns:\n")
print(colnames(seu_combined@meta.data))
cat("\nTreatment stage distribution:\n")
print(table(seu_combined@meta.data$treatment_stage))
cat("\nResponse distribution:\n")
print(table(seu_combined@meta.data$response))
cat("\nCell type distribution:\n")
print(table(seu_combined@meta.data$celltype))
seu_pre <- subset(seu_combined, subset = treatment_stage == "Pre-treatment")
seu_post <- subset(seu_combined, subset = treatment_stage == "Post-treatment")
seu_pre <- JoinLayers(seu_pre, assay = "RNA")
seu_post <- JoinLayers(seu_post, assay = "RNA")
saveRDS(seu_pre, file.path(output_dir, "seu_combined_harmony_annotated_Pre_treatment.rds"))
saveRDS(seu_post, file.path(output_dir, "seu_combined_harmony_annotated_Post_treatment.rds"))
cat("Data split finished, Pre-treatment cells:", ncol(seu_pre), "Post-treatment cells:", ncol(seu_post), "\n")
seu_epi_pre <- subset(seu_pre, subset = celltype == "Epithelial_cells")
seu_immune_pre <- subset(seu_pre, subset = celltype == "T_cells")
seu_epi_post <- subset(seu_post, subset = celltype == "Epithelial_cells")
seu_epi_pre <- JoinLayers(seu_epi_pre, assay = "RNA")
seu_immune_pre <- JoinLayers(seu_immune_pre, assay = "RNA")
seu_epi_post <- JoinLayers(seu_epi_post, assay = "RNA")
cat("\nLayer validation:\n")
cat("Epithelial pre layers:", Layers(seu_epi_pre, assay = "RNA"), "\n")
cat("T cell reference layers:", Layers(seu_immune_pre, assay = "RNA"), "\n")
cat("Post epithelial layers:", Layers(seu_epi_post, assay = "RNA"), "\n")
saveRDS(seu_epi_pre, file.path(output_dir, "seu_epi_pre.rds"))
saveRDS(seu_immune_pre, file.path(output_dir, "seu_immune_pre.rds"))
saveRDS(seu_epi_post, file.path(output_dir, "seu_epi_post.rds"))
cat("Cell subset extraction completed:\n")
cat("Pre-treatment epithelial cells:", ncol(seu_epi_pre), "\n")
cat("Pre-treatment T cell reference:", ncol(seu_immune_pre), "\n")
cat("Post-treatment epithelial cells:", ncol(seu_epi_post), "\n")
infercnv_dir <- file.path(output_dir, "InferCNV")
dir.create(infercnv_dir, showWarnings = FALSE)
if(ncol(seu_epi_pre) < 100) {
  warning("Low number of pre-treatment epithelial cells, may affect InferCNV results")
}
if(ncol(seu_immune_pre) < 20) {
  warning("Low number of T cell reference cells, may affect InferCNV results")
}
seu_epi_pre <- JoinLayers(seu_epi_pre, assay = "RNA")
counts_matrix <- LayerData(seu_epi_pre, assay = "RNA", layer = "counts")
seu_immune_pre <- JoinLayers(seu_immune_pre, assay = "RNA")
immune_counts <- LayerData(seu_immune_pre, assay = "RNA", layer = "counts")
common_genes <- intersect(rownames(counts_matrix), rownames(immune_counts))
combined_counts <- cbind(
  counts_matrix[common_genes, ],
  immune_counts[common_genes, ]
)
gencode <- read_tsv("hg38_gencode_v27.txt", col_names = c("gene", "chr", "start", "end"))
cat("First lines of gene position file:\n")
print(head(gencode))
genes_in_matrix <- rownames(combined_counts)
shared_genes <- intersect(genes_in_matrix, gencode$gene)
filtered_counts <- combined_counts[shared_genes, ]
original_colnames <- colnames(filtered_counts)
colnames(filtered_counts) <- original_colnames
cat("\nLast 10 column names of count matrix:\n")
if (ncol(filtered_counts) >= 10) {
  print(tail(colnames(filtered_counts), 10))
} else {
  print(colnames(filtered_counts))
}
write.table(filtered_counts,
            file.path(infercnv_dir, "infercnv_count_matrix.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
gencode_filtered <- gencode %>%
  filter(gene %in% shared_genes) %>%
  distinct(gene, .keep_all = TRUE)
write_tsv(gencode_filtered,
          file.path(infercnv_dir, "infercnv_gene_order.txt"),
          col_names = FALSE)
cat("\nFiltered count matrix dimension:", dim(filtered_counts), "\n")
cat("Filtered gene position file dimension:", dim(gencode_filtered), "\n")
cat("Number of shared genes:", length(shared_genes), "\n")
cell_annotations <- data.frame(
  cell_id = colnames(filtered_counts),
  cell_type = c(rep("Epithelial", ncol(seu_epi_pre)),
                rep("T_cells", ncol(seu_immune_pre)))
)
write.table(cell_annotations,
            file.path(infercnv_dir, "infercnv_cell_annotations.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("Cell annotations saved, cell type distribution:\n")
print(table(cell_annotations$cell_type))
cat("\nValidation between matrix columns and annotation cell_id:\n")
if (all(colnames(filtered_counts) == cell_annotations$cell_id)) {
  cat("Validation pass: column names match cell IDs perfectly\n")
} else {
  warning("Warning: column names mismatch with cell annotation IDs, please check!")
  mismatch_pos <- which(colnames(filtered_counts) != cell_annotations$cell_id)
  cat("Mismatch positions:", mismatch_pos, "\n")
}
count_matrix_genes <- nrow(filtered_counts)
cat("\nGene number in count matrix (infercnv_count_matrix.txt):", count_matrix_genes, "\n")
gene_order_genes <- nrow(gencode_filtered)
cat("Gene number in gene position file (infercnv_gene_order.txt):", gene_order_genes, "\n")
if (count_matrix_genes == gene_order_genes) {
  cat("Validation pass: gene counts consistent between matrix and gene position file\n")
} else {
  warning("Warning: gene number mismatch between matrix and gene position file, check filtering steps!")
}
