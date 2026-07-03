# CCND1 pCR Response Analysis Code

This repository contains the organized R scripts and lightweight metadata used for a CCND1-centered pCR versus NonpCR multi-omics workflow.

## Repository Layout

- `scripts/`: analysis scripts reorganized from `新整理代码.zip`.
- `metadata/`: lightweight sample/group metadata, gene-list metadata, and a checksum manifest for files in the source zip.

Large matrices, Seurat objects, and intermediate result tables are intentionally not committed to GitHub. See `metadata/data_manifest.tsv` for file sizes, MD5 checksums, and upload decisions.

## Scripts

| File | Purpose |
| --- | --- |
| `scripts/01_limma_DEA.R` | GSE25066 limma differential expression analysis. |
| `scripts/02_volcano_plot.R` | Volcano plot for differential expression results. |
| `scripts/03_DEG_heatmap.R` | DEG heatmap visualization. |
| `scripts/04_venn_intersection.R` | Venn/intersection step. The source zip file is empty and is retained as a placeholder. |
| `scripts/05_lasso_risk_model.R` | LASSO risk model training, risk score export, ROC, and coefficient plots. |
| `scripts/06_kegg_enrichment.R` | KEGG enrichment step. The source zip file is empty and is retained as a placeholder. |
| `scripts/07_gsea_hallmark_analysis.R` | MSigDB Hallmark GSEA analysis and pathway plots. |
| `scripts/08_scrna_merge_seurat.R` | Merge Seurat objects. |
| `scripts/09_scrna_qc_harmony_cluster.R` | Seurat QC, Harmony integration, clustering, and UMAP. |
| `scripts/10_scrna_celltype_annotation.R` | Marker-based scRNA-seq cell type annotation. |
| `scripts/11_scrna_infercnv_prepare.R` | Prepare inferCNV input files. |
| `scripts/12_scrna_infercnv_analysis_function_communication.R` | inferCNV analysis, malignant epithelial classification, and communication-related downstream functions. |

## Metadata

| File | Description |
| --- | --- |
| `metadata/GSE25066_group4.txt` | GSE25066 sample grouping metadata from the source zip. |
| `metadata/intersection_total.txt` | Gene/intersection metadata from the source zip. |
| `metadata/Diff-genes1.csv` | Differential gene list used by downstream enrichment/GSEA steps. |
| `metadata/data_manifest.tsv` | Full manifest of files in the source zip, including MD5 checksums and upload decisions. |

## Data Availability Notes

The following files were not uploaded because they are large data objects, analysis intermediates, or externally licensed reference resources:

- `1/GSE25066_symbol_avg_filtered_ordered.csv` (about 109 MB)
- `6/seu_combined.rds` (about 389 MB)
- intermediate DEG/LASSO/inferCNV input tables listed in `metadata/data_manifest.tsv`
- `5/h.all.v2024.1.Hs.symbols.gmt`, which should be downloaded from MSigDB rather than redistributed here

For manuscript or review use, place the required input files next to the relevant script or edit the file paths. Public datasets should be cited by accession and download date; generated or private data should be deposited in an appropriate data repository rather than committed directly to GitHub.

## Reproducibility Notes

Before claiming end-to-end reproducibility, record:

- R version and Bioconductor version
- package versions
- dataset accessions and download dates
- random seed
- input file checksums
- output file checksums

The scripts are preserved as analysis code from the organized source zip. Some steps still require local input files and may need path updates before execution.
