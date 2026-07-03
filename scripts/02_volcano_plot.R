library(ggplot2)
library(ggrepel)
library(tidyverse)
library(data.table)
data <- fread('GSE25066_diff.txt')
log2FCfilter = 0.585
log2FCcolor = 0.585
colnames(data)[2] <- 'log2FC'
colnames(data)[1] <- 'gene'
index = data$adj.P.Val <0.05 & abs(data$log2FC) > log2FCfilter
data$group <- 0
data$group[index & data$log2FC>0] = 1
data$group[index & data$log2FC<0] = -1
data$group <- factor(data$group,levels = c(1,0,-1),labels =c("UP","NOT","DOWN") )
data <- as.data.frame(data)
rownames(data) <- data$gene
data <- data[, -which(colnames(data) == "gene"), drop = FALSE]
write.table(data,file="GSE25066_diff1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
data$gene <- rownames(data)
p <- ggplot(data = data, aes(x = log2FC, y = -log10(adj.P.Val), color = group)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values =  c("#E64B35FF", "grey50", "#4DBBD5FF"), labels = c("UP ", "NOT", "DOWN ")) +
  labs(x = "log2 (fold change)", y = "-log10 (adj.P.Val)") +
  theme(plot.title = element_text(hjust = 0.4)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = c(-log2FCfilter, log2FCfilter), lty = 4, lwd = 0.6, alpha = 0.8) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "top")
pdf("volcano_plot.pdf", width = 8, height = 8)
print(p)
dev.off()
tiff("volcano_plot.tiff", width = 8, height = 8, units = "in", res = 300, compression = "lzw")
print(p)
dev.off()
