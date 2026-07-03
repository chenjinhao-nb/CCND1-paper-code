library(Seurat)
library(ggplot2)
library(patchwork)
seu_combined <- readRDS("seu_combined_harmony.rds")
breast_cell_markers <- list(
  "Epithelial_cells" = c("EPCAM", "KRT8", "KRT18"),
  "Endothelial_cells" = c("VWF", "PLVAP", "PECAM1"),
  "Pericytes" = c("RGS5", "ACTA2", "MCAM"),
  "Fibroblasts" = c("COL1A1", "DCN", "LUM"),
  "Myeloid_cells" = c("LYZ", "CD68", "CD14"),
  "T_cells" = c("CD3D", "CD3E", "CD2"),
  "B_cells" = c("MS4A1", "CD79A", "CD79B"),
  "Plasma_cells" = c("IGHG1", "JCHAIN", "IGKC"),
  "NK_cells" = c("KLRD1", "GNLY", "NKG7"),
  "Mast_cells" = c("TPSAB1", "TPSB2", "KIT"),
  "Proliferating_cells" = c("MKI67", "PCNA", "TOP2A")
)
p_dot <- DotPlot(seu_combined, features = unique(unlist(breast_cell_markers)), assay = "RNA", cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Gene Expression Across Clusters")
ggsave("marker_gene_dotplot.pdf", p_dot, width = 12, height = 8)
expression_data <- DotPlot(seu_combined, features = unique(unlist(breast_cell_markers)), assay = "RNA", cluster.idents = TRUE)$data
cluster_ids <- levels(expression_data$id)
celltype_scores <- matrix(0, nrow = length(cluster_ids), ncol = length(breast_cell_markers),
                          dimnames = list(cluster_ids, names(breast_cell_markers)))
for (celltype in names(breast_cell_markers)) {
  markers <- breast_cell_markers[[celltype]]
  for (cluster in cluster_ids) {
    cluster_data <- expression_data[expression_data$id == cluster & expression_data$features.plot %in% markers, ]
    score <- mean(cluster_data$avg.exp.scaled * cluster_data$pct.exp / 100, na.rm = TRUE)
    celltype_scores[cluster, celltype] <- ifelse(is.na(score), 0, score)
  }
}
cluster_annotations <- data.frame(
  Cluster = cluster_ids,
  Predicted_CellType = colnames(celltype_scores)[apply(celltype_scores, 1, which.max)],
  Max_Score = apply(celltype_scores, 1, max)
)
score_df <- as.data.frame(celltype_scores)
colnames(score_df) <- paste0("Score_", colnames(score_df))
cluster_annotations <- cbind(cluster_annotations, score_df)
write.csv(cluster_annotations, "automated_celltype_annotations.csv", row.names = FALSE)
cluster_celltype_manual <- c(
  "10" = "B_cells",
  "16" = "Mast_cells",
  "17" = "Epithelial_cells",
  "18" = "Plasma_cells",
  "11" = "Epithelial_cells",
  "8" = "Mast_cells",
  "5" = "Epithelial_cells",
  "1" = "Epithelial_cells",
  "4" = "Epithelial_cells",
  "14" = "Epithelial_cells",
  "2" = "Myeloid_cells",
  "19" = "Myeloid_cells",
  "21" = "Plasma_cells",
  "7" = "T_cells",
  "0" = "T_cells",
  "3" = "T_cells",
  "13" = "NK_cells",
  "15" = "Proliferating_cells",
  "20" = "Proliferating_cells",
  "9" = "Endothelial_cells",
  "6" = "Fibroblasts",
  "12" = "Pericytes"
)
print("Manual annotation results:")
print(cluster_celltype_manual)
seu_combined[['celltype']] <- unname(cluster_celltype_manual[as.character(seu_combined@meta.data$seurat_clusters)])
p_umap <- DimPlot(seu_combined, reduction = 'umap', group.by = 'celltype', label = TRUE, label.size = 3, pt.size = 0.1, repel = TRUE) +
  ggtitle("UMAP by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("umap_by_celltype.pdf", p_umap, width = 10, height = 8)
celltype_counts <- table(seu_combined$celltype)
print("Cell type distribution:")
print(celltype_counts)
write.csv(celltype_counts, "celltype_counts.csv", row.names = FALSE)
celltype_df <- as.data.frame(celltype_counts)
names(celltype_df) <- c("CellType", "Count")
celltype_df$Percentage <- round(celltype_df$Count / sum(celltype_df$Count) * 100, 1)
p_pie <- ggplot(celltype_df, aes(x = "", y = Count, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Cell Type Composition")
ggsave("celltype_composition.pdf", p_pie, width = 8, height = 6)
p_dot_celltype <- DotPlot(seu_combined, features = unique(unlist(breast_cell_markers)), assay = "RNA", group.by = 'celltype') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Gene Expression by Cell Type")
ggsave("marker_gene_dotplot_by_celltype.pdf", p_dot_celltype, width = 12, height = 8)
print(p_umap)
print(p_dot_celltype)
saveRDS(seu_combined, file = "seu_combined_harmony_annotated.rds")
