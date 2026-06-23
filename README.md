# CCND1 pCR Response Analysis Code

This repository records the R analysis scripts used for a CCND1-centered pCR versus NonpCR multi-omics workflow.

## Contents

- `scripts/R1.txt`: GEO microarray preprocessing for GSE25066, normalization, annotation, PCA, limma differential expression, and volcano plot.
- `scripts/02_Functional_EnrichmentGO_enrichment.txt`: GO enrichment and visualization.
- `scripts/02_Functional_EnrichmentKEGG_enrichment.txt`: KEGG enrichment, pathway network, and pathview outputs.
- `scripts/02_Functional_EnrichmentHallmark_GSEA.txt`: MSigDB Hallmark GSEA and pathway plots.
- `scripts/03_DEG_VisualizationDEG_Heatmap.txt`: DEG heatmap generation.
- `scripts/04_LASSO_Signature_ModelLASSO_data_prep.txt`: LASSO input matrix preparation.
- `scripts/04_LASSO_Signature_ModelLASSO_model_train.txt`: LASSO model training, risk score, ROC, and coefficient export.
- `scripts/05_Serum_MetabolomicsMetabolomics_OPLSDA.txt`: Serum metabolomics OPLS-DA and MetaboAnalystR analysis.
- `scripts/06_scRNAseq_AnalysisscRNA_QC_Harmony.txt`: Seurat QC, Harmony integration, clustering, and UMAP.
- `scripts/06_scRNAseq_AnalysisscRNA_Cell_Annotation.txt`: Marker-based cell type annotation.
- `scripts/06_scRNAseq_AnalysisscRNA_inferCNV.txt`: inferCNV-based malignant epithelial cell classification.
- `scripts/06_scRNAseq_AnalysisCCND1_Correlation_Analysis.txt`: CCND1 expression and correlation analysis in malignant epithelial cells.
- `archive/original-code.zip`: Original code archive supplied before repository organization.

## Reproducibility Notes

The scripts expect local input files such as expression matrices, sample grouping files, RDS objects, annotation tables, and MSigDB GMT files. Before running, place the required input files next to the relevant script or edit the file paths.

Some extracted scripts appear to contain typos or truncation. Review syntax and object names before claiming end-to-end reproducibility.

Recommended metadata to record when running these scripts:

- R version and Bioconductor version
- package versions
- source dataset accession and download date
- random seed
- input file checksums
- generated output checksums
