library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(gridExtra)
if (!dir.exists("hall")) {
  dir.create("hall")
}
diff_data <- read.delim("GSE25066_diff.txt", header = TRUE, stringsAsFactors = FALSE)
if(any(duplicated(diff_data$gene_symbol))) {
  diff_data <- diff_data[order(diff_data$gene_symbol, -abs(diff_data$logFC)), ]
  diff_data <- diff_data[!duplicated(diff_data$gene_symbol), ]
}
gene_list <- diff_data$logFC
names(gene_list) <- diff_data$gene_symbol
gene_list <- sort(gene_list, decreasing = TRUE)
set.seed(123)
hallmark_gene_sets <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
gsea_res <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = hallmark_gene_sets,
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE
)
gsea_res@result$Description <- gsub("HALLMARK_", "", gsea_res@result$Description)
write.csv(gsea_res@result, "hall/GSEA_Hallmark_results.csv", row.names = FALSE)
sig_res <- subset(gsea_res@result, p.adjust < 0.05)
if(nrow(sig_res) > 0) {
  dotplot_data <- head(sig_res[order(-abs(sig_res$NES)), ], 10)
  p_dot <- ggplot(dotplot_data,
                  aes(x = NES, y = reorder(Description, NES),
                      color = -log10(p.adjust), size = setSize)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(adj.p)") +
    labs(x = "Normalized Enrichment Score (NES)",
         y = NULL,
         title = "Top 10 Enriched Hallmark Pathways",
         size = "Gene Set Size") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))
  pdf("hall/GSEA_Hallmark_dotplot.pdf", width = 10, height = 6)
  print(p_dot)
  dev.off()
  tiff("hall/GSEA_Hallmark_dotplot.tiff", width = 10, height = 6, units = "in", res = 300)
  print(p_dot)
  dev.off()
  if(nrow(sig_res) >= 1) {
    top_pathway <- sig_res$ID[1]
    p_enrich <- enrichplot::gseaplot2(
      gsea_res,
      geneSetID = top_pathway,
      title = sig_res$Description[1],
      pvalue_table = TRUE,
      color = "firebrick",
      base_size = 12
    )
    pdf("hall/Top_Hallmark_pathway.pdf", width = 14, height = 6)
    print(p_enrich)
    dev.off()
    tiff("hall/Top_Hallmark_pathway.tiff", width = 14, height = 6, units = "in", res = 300)
    print(p_enrich)
    dev.off()
  }
  if(nrow(sig_res) >= 3) {
    top3 <- sig_res$ID[1:3]
    p_multi <- enrichplot::gseaplot2(
      gsea_res,
      geneSetID = top3,
      pvalue_table = FALSE,
      base_size = 10,
      color = c("#E41A1C", "#377EB8", "#4DAF4A")
    )
    pdf("hall/Top3_Hallmark_pathways.pdf", width = 14, height = 8)
    print(p_multi)
    dev.off()
    tiff("hall/Top3_Hallmark_pathways.tiff", width = 14, height = 8, units = "in", res = 300)
    print(p_multi)
    dev.off()
  }
  if(nrow(sig_res) >= 1) {
    p_net <- cnetplot(gsea_res,
                      categorySize = "pvalue",
                      foldChange = gene_list,
                      showCategory = min(5, nrow(sig_res)),
                      colorEdge = TRUE,
                      node_label = "all") +
      ggtitle("Gene-Concept Network (Hallmark Pathways)") +
      theme(legend.position = "right")
    pdf("hall/Pathway_network_Hallmark.pdf", width = 12, height = 10)
    print(p_net)
    dev.off()
    tiff("hall/Pathway_network_Hallmark.tiff", width = 12, height = 10, units = "in", res = 300)
    print(p_net)
    dev.off()
  }
} else {
  message("No significantly enriched Hallmark pathways (adj.p < 0.05)")
}
if(nrow(sig_res) >= 21) {
  selected <- sig_res$ID[c(3, 6, 21)]
  p_multi <- enrichplot::gseaplot2(
    gsea_res,
    geneSetID = selected,
    pvalue_table = FALSE,
    base_size = 10,
    color = c("#E41A1C", "#377EB8", "#4DAF4A")
  )
  pdf("3_pathways.pdf", width = 14, height = 8)
  print(p_multi)
  dev.off()
  tiff("3_pathways.tiff", width = 14, height = 8, units = "in", res = 300)
  print(p_multi)
  dev.off()
} else {
  message("Less than 21 significant pathways, cannot generate target plot")
}
if(nrow(sig_res) >= 22) {
  top_pathway <- sig_res$ID[22]
  p_enrich <- enrichplot::gseaplot2(
    gsea_res,
    geneSetID = top_pathway,
    title = sig_res$Description[22],
    pvalue_table = FALSE,
    color = "firebrick",
    base_size = 12
  )
  pdf("hall/22_pathway.pdf", width = 14, height = 6)
  print(p_enrich)
  dev.off()
  tiff("hall/22_pathway.tiff", width = 14, height = 6, units = "in", res = 300)
  print(p_enrich)
  dev.off()
}
