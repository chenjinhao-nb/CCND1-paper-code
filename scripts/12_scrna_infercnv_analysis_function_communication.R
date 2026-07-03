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
library(tibble)
library(pheatmap)
seu_files <- c(
  "CCND1_Analysis_Output/seu_immune_pre.rds",
  "CCND1_Analysis_Output/seu_epi_pre.rds",
  "CCND1_Analysis_Output/seu_epi_post.rds",
  "CCND1_Analysis_Output/seu_combined_harmony_annotated_Pre_treatment.rds",
  "CCND1_Analysis_Output/seu_combined_harmony_annotated_Post_treatment.rds"
)
seu_immune_pre <- readRDS(seu_files[1])
seu_epi_pre <- readRDS(seu_files[2])
seu_epi_post <- readRDS(seu_files[3])
seu_combined_pre <- readRDS(seu_files[4])
seu_combined_post <- readRDS(seu_files[5])
infercnv_gene_order <- read.table("CCND1_Analysis_Output/InferCNV/infercnv_gene_order.txt", header = FALSE, sep = "\t")
infercnv_count_matrix <- read.table("CCND1_Analysis_Output/InferCNV/infercnv_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
infercnv_cell_annotations <- read.table("CCND1_Analysis_Output/InferCNV/infercnv_cell_annotations.txt", header = FALSE, sep = "\t")
colnames(infercnv_cell_annotations) <- c("cell_id", "cell_type")
output_dir <- "InferCNV_Results"
infercnv_dir <- file.path(output_dir, "InferCNV_Analysis")
dir.create(infercnv_dir, recursive = TRUE, showWarnings = FALSE)
temp_annot_file <- file.path(infercnv_dir, "temp_cell_annotations.txt")
write.table(infercnv_cell_annotations, temp_annot_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "CCND1_Analysis_Output/InferCNV/infercnv_count_matrix.txt",
  annotations_file = temp_annot_file,
  delim = "\t",
  gene_order_file = "CCND1_Analysis_Output/InferCNV/infercnv_gene_order.txt",
  ref_group_names = c("T_cells")
)
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = infercnv_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = 6,
  analysis_mode = "subclusters"
)
if(file.exists(temp_annot_file)) file.remove(temp_annot_file)
cat("Start saving analysis results...\n")
saveRDS(infercnv_obj, file = file.path(infercnv_dir, "infercnv_final_object.rds"))
cat("- Full infercnv object saved: infercnv_final_object.rds\n")
if (!is.null(infercnv_obj@expr.data)) {
  write.table(infercnv_obj@expr.data,
              file = file.path(infercnv_dir, "infercnv.observations.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
  cat("- Observation matrix saved: infercnv.observations.txt\n")
  observations_file <- file.path(infercnv_dir, "infercnv.observations.txt")
  if (file.exists(observations_file)) {
    file.copy(observations_file,
              file.path(infercnv_dir, "infercnv_expr_matrix.txt"),
              overwrite = TRUE)
    cat("- Visualization matrix saved: infercnv_expr_matrix.txt\n")
  }
}
write.table(infercnv_obj@gene_order,
            file = file.path(infercnv_dir, "infercnv_gene_positions.txt"),
            sep = "\t", quote = FALSE)
cat("- Gene position file saved: infercnv_gene_positions.txt\n")
obs_file <- file.path(infercnv_dir, "infercnv.observations.txt")
if(file.exists(obs_file)) {
  cnv_scores <- read.table(obs_file, header = TRUE, row.names = 1)
  cell_cnv_means <- colMeans(abs(cnv_scores))
  epithelial_cells <- infercnv_cell_annotations$cell_id[infercnv_cell_annotations$cell_type == "Epithelial"]
  epithelial_cnv <- cell_cnv_means[names(cell_cnv_means) %in% epithelial_cells]
  cnv_threshold <- quantile(epithelial_cnv, 0.75)
  cat("CNV threshold:", cnv_threshold, "\n")
  malignant_cells <- names(epithelial_cnv[epithelial_cnv > cnv_threshold])
  seu_epi_pre$Cell_Status <- ifelse(colnames(seu_epi_pre) %in% malignant_cells, "Malignant", "Non_Malignant")
  cnv_plot_data <- data.frame(
    cell_id = names(cell_cnv_means),
    cnv_score = cell_cnv_means,
    cell_type = infercnv_cell_annotations$cell_type[match(names(cell_cnv_means), infercnv_cell_annotations$cell_id)]
  )
  p_cnv <- ggplot(cnv_plot_data, aes(x = cell_type, y = cnv_score, fill = cell_type)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.7) +
    geom_hline(yintercept = cnv_threshold, linetype = "dashed", color = "red") +
    labs(title = "CNV Scores by Cell Type",
         subtitle = paste("Threshold:", round(cnv_threshold, 3)),
         x = "Cell Type", y = "CNV Score") +
    theme_minimal()
  ggsave(file.path(infercnv_dir, "CNV_score_distribution.pdf"), p_cnv, width = 8, height = 6)
  cnv_results <- data.frame(
    cell_id = names(cell_cnv_means),
    cnv_score = cell_cnv_means,
    cell_type = cnv_plot_data$cell_type,
    is_malignant = names(cell_cnv_means) %in% malignant_cells
  )
  write.csv(cnv_results, file.path(output_dir, "cnv_scores_detailed.csv"))
  cat("Malignant cell identification finished:\n")
  cat("Malignant cells:", sum(seu_epi_pre$Cell_Status == "Malignant"), "\n")
  cat("Non-malignant epithelial cells:", sum(seu_epi_pre$Cell_Status == "Non_Malignant"), "\n")
} else {
  warning("InferCNV output missing, skip CNV scoring analysis")
  seu_epi_pre$Cell_Status <- "Unknown"
}
saveRDS(seu_epi_pre, file.path(output_dir, "seu_epi_pre_with_CNV.rds"))
cat("Start STC2+ tumor cell functional analysis...\n")
if("STC2" %in% rownames(seu_epi_pre)) {
  stc2_expression <- GetAssayData(seu_epi_pre, slot = "data")["STC2", ]
  stc2_threshold <- quantile(stc2_expression, 0.75)
  seu_epi_pre$STC2_Status <- ifelse(stc2_expression > stc2_threshold, "STC2_High", "STC2_Low")
  seu_epi_pre$Analysis_Group <- case_when(
    seu_epi_pre$Cell_Status == "Malignant" & seu_epi_pre$STC2_Status == "STC2_High" ~ "STC2_Malignant",
    seu_epi_pre$Cell_Status == "Non_Malignant" ~ "Normal_Epi",
    TRUE ~ "Other_Malignant"
  )
  if("response" %in% colnames(seu_combined_pre@meta.data)) {
    seu_epi_pre$response <- seu_combined_pre$response[match(colnames(seu_epi_pre), colnames(seu_combined_pre))]
  } else {
    warning("Response metadata not found, skip response-based subgroup analysis")
  }
  glycolysis_genes <- c("HK2", "GP1", "ALDOA", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA")
  glycolysis_genes <- glycolysis_genes[glycolysis_genes %in% rownames(seu_epi_pre)]
  cell_cycle_genes <- c("MKI67", "PCNA", "TOP2A", "BIRC5")
  cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes %in% rownames(seu_epi_pre)]
  drug_resistance_genes <- c("ABCB1", "ABCG2", "ABCC1")
  drug_resistance_genes <- drug_resistance_genes[drug_resistance_genes %in% rownames(seu_epi_pre)]
  if(length(glycolysis_genes) > 3) {
    seu_epi_pre <- AddModuleScore(seu_epi_pre, features = list(glycolysis_genes), name = "Glycolysis_Score")
  }
  if(length(cell_cycle_genes) > 3) {
    seu_epi_pre <- AddModuleScore(seu_epi_pre, features = list(cell_cycle_genes), name = "CellCycle_Score")
  }
  if(length(drug_resistance_genes) > 1) {
    seu_epi_pre <- AddModuleScore(seu_epi_pre, features = list(drug_resistance_genes), name = "DrugResistance_Score")
  }
  if("response" %in% colnames(seu_epi_pre@meta.data)) {
    score_columns <- c("Analysis_Group", "response")
    if("Glycolysis_Score1" %in% colnames(seu_epi_pre@meta.data)) score_columns <- c(score_columns, "Glycolysis_Score1")
    if("CellCycle_Score1" %in% colnames(seu_epi_pre@meta.data)) score_columns <- c(score_columns, "CellCycle_Score1")
    if("DrugResistance_Score1" %in% colnames(seu_epi_pre@meta.data)) score_columns <- c(score_columns, "DrugResistance_Score1")
    functional_data <- FetchData(seu_epi_pre, vars = score_columns)
    if(ncol(functional_data) > 2) {
      functional_long <- functional_data %>%
        pivot_longer(cols = -c(Analysis_Group, response), names_to = "Pathway", values_to = "Score")
      p_functional_response <- ggplot(functional_long, aes(x = Analysis_Group, y = Score, fill = response)) +
        geom_boxplot(alpha = 0.7) +
        facet_wrap(~ Pathway, scales = "free_y", nrow = 1) +
        labs(title = "Functional Pathway Scores by Cell Group and Response",
             x = "Cell Group", y = "Module Score") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(output_dir, "Functional_scores_by_response.pdf"), p_functional_response, width = 12, height = 6)
    }
  }
  Idents(seu_epi_pre) <- "Analysis_Group"
  if(exists("seu_epi_pre$response") &&
     sum(seu_epi_pre$Analysis_Group == "STC2_Malignant" & seu_epi_pre$response == "pCR") > 5 &&
     sum(seu_epi_pre$Analysis_Group == "Normal_Epi" & seu_epi_pre$response == "pCR") > 5) {
    de_genes_pCR <- FindMarkers(
      subset(seu_epi_pre, subset = response == "pCR"),
      ident.1 = "STC2_Malignant",
      ident.2 = "Normal_Epi",
      min.pct = 0.25,
      logfc.threshold = 0.5
    )
    write.csv(de_genes_pCR, file.path(output_dir, "STC2_Malignant_vs_Normal_DE_genes_pCR.csv"))
  }
  if(exists("seu_epi_pre$response") &&
     sum(seu_epi_pre$Analysis_Group == "STC2_Malignant" & seu_epi_pre$response == "NonpCR") > 5 &&
     sum(seu_epi_pre$Analysis_Group == "Normal_Epi" & seu_epi_pre$response == "NonpCR") > 5) {
    de_genes_NonpCR <- FindMarkers(
      subset(seu_epi_pre, subset = response == "NonpCR"),
      ident.1 = "STC2_Malignant",
      ident.2 = "Normal_Epi",
      min.pct = 0.25,
      logfc.threshold = 0.5
    )
    write.csv(de_genes_NonpCR, file.path(output_dir, "STC2_Malignant_vs_Normal_DE_genes_NonpCR.csv"))
  }
  cat("Functional analysis completed, outputs stored in:", output_dir, "\n")
} else {
  warning("STC2 gene not detected in expression matrix, skip STC2-related analysis")
}
cat("Preparing input files for cell-cell communication analysis...\n")
seu_pre_updated <- seu_combined_pre
seu_pre_updated$Cell_Status <- ifelse(colnames(seu_pre_updated) %in% colnames(seu_epi_pre), seu_epi_pre$Cell_Status, "Other")
comm_dir <- file.path(output_dir, "CellCommunication")
dir.create(comm_dir, showWarnings = FALSE)
writeMM(GetAssayData(seu_pre_updated, slot = "counts"), file.path(comm_dir, "counts_matrix.mtx"))
write.table(
  data.frame(gene = rownames(seu_pre_updated)),
  file.path(comm_dir, "gene_names.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
if("celltype" %in% colnames(seu_pre_updated@meta.data)) {
  comm_annot <- data.frame(
    Cell = colnames(seu_pre_updated),
    CellType = seu_pre_updated$celltype
  )
} else {
  comm_annot <- data.frame(
    Cell = colnames(seu_pre_updated),
    CellType = "Unknown"
  )
}
write.csv(comm_annot, file.path(comm_dir, "cell_annotations.csv"), row.names = FALSE)
if("response" %in% colnames(seu_pre_updated@meta.data)) {
  response_annot <- data.frame(
    Cell = colnames(seu_pre_updated),
    Response = seu_pre_updated$response,
    Sample = seu_pre_updated$sample_id
  )
} else {
  response_annot <- data.frame(
    Cell = colnames(seu_pre_updated),
    Response = "Unknown",
    Sample = "Unknown"
  )
}
write.csv(response_annot, file.path(comm_dir, "response_info.csv"), row.names = FALSE)
cat("CellPhoneDB input files generated:\n")
cat("Count matrix:", file.path(comm_dir, "counts_matrix.mtx"), "\n")
cat("Cell annotation:", file.path(comm_dir, "cell_annotations.csv"), "\n")
cat("Response metadata:", file.path(comm_dir, "response_info.csv"), "\n")
cat("Start time trajectory & dynamic change analysis...\n")
if("STC2" %in% rownames(seu_epi_post) && exists("stc2_threshold")) {
  stc2_expression_post <- GetAssayData(seu_epi_post, slot = "data")["STC2", ]
  seu_epi_post$STC2_Status <- ifelse(stc2_expression_post > stc2_threshold, "STC2_High", "STC2_Low")
  if("response" %in% colnames(seu_combined_post@meta.data)) {
    seu_epi_post$response <- seu_combined_post$response[match(colnames(seu_epi_post), colnames(seu_combined_post))]
  }
  if(exists("glycolysis_genes") && length(glycolysis_genes) > 3) {
    seu_epi_post <- AddModuleScore(seu_epi_post, features = list(glycolysis_genes), name = "Glycolysis_Score")
  }
  if(exists("cell_cycle_genes") && length(cell_cycle_genes) > 3) {
    seu_epi_post <- AddModuleScore(seu_epi_post, features = list(cell_cycle_genes), name = "CellCycle_Score")
  }
  if(exists("seu_epi_pre$response") && exists("seu_epi_post$response")) {
    pre_composition <- table(seu_epi_pre$STC2_Status, seu_epi_pre$response)
    pre_composition_prop <- prop.table(pre_composition, margin = 2)
    post_composition <- table(seu_epi_post$STC2_Status, seu_epi_post$response)
    post_composition_prop <- prop.table(post_composition, margin = 2)
    comp_data <- rbind(
      as.data.frame(pre_composition_prop) %>%
        mutate(Timepoint = "Pre", Proportion = Freq) %>%
        select(Timepoint, STC2_Status = Var1, Response = Var2, Proportion),
      as.data.frame(post_composition_prop) %>%
        mutate(Timepoint = "Post", Proportion = Freq) %>%
        select(Timepoint, STC2_Status = Var1, Response = Var2, Proportion)
    )
    p_composition_response <- ggplot(comp_data, aes(x = Timepoint, y = Proportion, fill = STC2_Status)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ Response) +
      labs(title = "STC2+ Cell Proportion Changes After NAC by Response",
           x = "Timepoint", y = "Proportion") +
      theme_minimal()
    ggsave(file.path(output_dir, "STC2_cell_proportion_dynamics_by_response.pdf"), p_composition_response, width = 10, height = 6)
  }
  score_columns_dynamics <- c()
  if("Glycolysis_Score1" %in% colnames(seu_epi_pre@meta.data)) score_columns_dynamics <- c(score_columns_dynamics, "Glycolysis_Score1")
  if("CellCycle_Score1" %in% colnames(seu_epi_pre@meta.data)) score_columns_dynamics <- c(score_columns_dynamics, "CellCycle_Score1")
  if(length(score_columns_dynamics) > 0 && exists("seu_epi_pre$response") && exists("seu_epi_post$response")) {
    functional_dynamics_pre <- FetchData(seu_epi_pre, vars = score_columns_dynamics) %>%
      mutate(Timepoint = "Pre", Response = seu_epi_pre$response)
    functional_dynamics_post <- FetchData(seu_epi_post, vars = score_columns_dynamics) %>%
      mutate(Timepoint = "Post", Response = seu_epi_post$response)
    functional_dynamics <- rbind(functional_dynamics_pre, functional_dynamics_post)
    functional_dynamics_long <- functional_dynamics %>%
      pivot_longer(cols = all_of(score_columns_dynamics), names_to = "Pathway", values_to = "Score")
    p_dynamics_response <- ggplot(functional_dynamics_long, aes(x = Timepoint, y = Score, fill = Timepoint)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.2, alpha = 0.7) +
      facet_grid(Pathway ~ Response, scales = "free_y") +
      labs(title = "Functional Pathway Dynamics After NAC by Response",
           x = "Timepoint", y = "Module Score") +
      theme_minimal()
    ggsave(file.path(output_dir, "Functional_dynamics_by_response.pdf"), p_dynamics_response, width = 12, height = 8)
  }
} else {
  warning("STC2 missing in post-treatment dataset or threshold undefined, skip dynamic trajectory analysis")
}
cat("Saving final objects and generating analysis summary...\n")
saveRDS(seu_epi_pre, file.path(output_dir, "final_seu_epi_pre_analyzed.rds"))
saveRDS(seu_epi_post, file.path(output_dir, "final_seu_epi_post_analyzed.rds"))
cat("\n=== Analysis Summary Report ===\n")
cat("Output directory:", output_dir, "\n")
cat("Total pre-treatment epithelial cells:", ncol(seu_epi_pre), "\n")
if(exists("malignant_cells")) {
  cat("Identified malignant cells:", sum(seu_epi_pre$Cell_Status == "Malignant"), "\n")
  cat("Non-malignant epithelial cells:", sum(seu_epi_pre$Cell_Status == "Non_Malignant"), "\n")
}
if("STC2_Status" %in% colnames(seu_epi_pre@meta.data)) {
  cat("STC2 high-expression cells:", sum(seu_epi_pre$STC2_Status == "STC2_High"), "\n")
}
cat("Total post-treatment epithelial cells:", ncol(seu_epi_post), "\n")
if("response" %in% colnames(seu_epi_pre@meta.data)) {
  cat("\nPre-treatment cell count grouped by response & cell status:\n")
  print(table(seu_epi_pre$response, seu_epi_pre$Cell_Status))
}
cat("Analysis finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nKey output paths:\n")
cat("- InferCNV CNV results:", infercnv_dir, "\n")
cat("- Functional analysis figures:", file.path(output_dir, "*.pdf"), "\n")
cat("- CellPhoneDB input files:", comm_dir, "\n")
cat("- Final analyzed Seurat objects:", file.path(output_dir, "final_seu_epi_*_analyzed.rds"), "\n")
