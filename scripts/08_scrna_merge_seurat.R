library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(AUCell)
library(msigdbr)
data_dir <- "./scRNA_input/"
sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
sample_names <- basename(sample_dirs)
sample_metadata <- data.frame(
  sample_id = sample_names,
  treatment_stage = ifelse(grepl("GSM621332[2-9]", sample_names),
                           "Pre-treatment", "Post-treatment"),
  response = ifelse(grepl("GSM621332[2-4]|GSM621333[0-1]", sample_names),
                    "pCR", "NonpCR"),
  stringsAsFactors = FALSE
)
seu_list <- list()
for (i in seq_along(sample_names)) {
  cat("Processing sample:", sample_names[i], "\n")
  counts_data <- Read10X(data.dir = sample_dirs[i])
  seu_obj <- CreateSeuratObject(
    counts = counts_data,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200
  )
  meta_row <- sample_metadata[sample_metadata$sample_id == sample_names[i], ]
  seu_obj$treatment_stage <- meta_row$treatment_stage
  seu_obj$response <- meta_row$response
  seu_obj$sample_id <- sample_names[i]
  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
  seu_obj[["percent.rb"]] <- PercentageFeatureSet(seu_obj, pattern = "^RP[SL]")
  seu_list[[i]] <- seu_obj
}
names(seu_list) <- sample_names
seu_combined <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = sample_names
)
saveRDS(seu_combined, file = "seu_combined.rds")
cell_counts <- seu_combined@meta.data %>%
  group_by(sample_id) %>%
  summarise(cell_number = n(), .groups = "drop") %>%
  arrange(sample_id)
print("Cell count per sample:")
print(cell_counts)
write.csv(cell_counts,
          file = "sample_cell_counts.csv",
          row.names = FALSE)
