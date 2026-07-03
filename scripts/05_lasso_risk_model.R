if (!require("glmnet")) install.packages("glmnet")
if (!require("pROC")) install.packages("pROC")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(glmnet)
library(pROC)
library(RColorBrewer)
data <- read.csv("GSE25066-LASSO1.csv", header = TRUE, check.names = FALSE)
y <- as.numeric(trimws(data[, 2]) == "NonpCR")
if (length(unique(y)) < 2) {
  stop("Response variable must contain two groups (pCR and NonpCR), please check input data.")
}
x <- as.matrix(data[, 3:ncol(data)])
if (any(is.na(x)) || any(is.na(y))) {
  stop("Missing values detected in dataset, please clean data before running.")
}
set.seed(123)
cv_fit <- cv.glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "auc"
)
risk_scores <- predict(cv_fit, newx = x, s = "lambda.min", type = "response")
risk_scores_vector <- as.vector(risk_scores)
selected_genes <- coef(cv_fit, s = "lambda.min")[which(coef(cv_fit, s = "lambda.min") != 0),]
selected_genes <- selected_genes[-1]
selected_gene_names <- names(selected_genes)
data_with_risk <- cbind(
  data[, 1:2],
  RiskScore = risk_scores_vector,
  data[, selected_gene_names, drop = FALSE]
)
write.csv(data_with_risk, "GSE25066-LASSO2.csv", row.names = FALSE)
palette <- brewer.pal(9, "Set1")
par(mar = c(5, 5, 4, 2) + 0.1,
    cex.lab = 1.2,
    cex.axis = 1,
    family = "sans")
pdf("CV_Curve_Optimized.pdf", width = 8, height = 8)
plot(cv_fit,
     main = "",
     col = palette[2],
     cex.main = 1.5,
     xlab = expression(log(lambda)),
     ylab = "AUC (10-fold CV)")
abline(v = log(cv_fit$lambda.min), col = palette[1], lty = 2, lwd = 2)
abline(v = log(cv_fit$lambda.1se), col = palette[3], lty = 3, lwd = 2)
legend("topright",
       legend = c("Optimal lambda (max AUC)", "1-SE rule lambda"),
       col = palette[c(1,3)], lty = c(2,3), lwd = 2,
       cex = 0.9, bty = "n")
dev.off()
tiff("CV_Curve_Optimized.tif", width = 8, height = 8, units = "in", res = 300, compression = "lzw")
plot(cv_fit,
     main = "",
     col = palette[2],
     cex.main = 1.5,
     xlab = expression(log(lambda)),
     ylab = "AUC (10-fold CV)")
abline(v = log(cv_fit$lambda.min), col = palette[1], lty = 2, lwd = 2)
abline(v = log(cv_fit$lambda.1se), col = palette[3], lty = 3, lwd = 2)
legend("topright",
       legend = c("Optimal lambda (max AUC)", "1-SE rule lambda"),
       col = palette[c(1,3)], lty = c(2,3), lwd = 2,
       cex = 0.9, bty = "n")
dev.off()
pdf("Coefficient_Path_Optimized.pdf", width = 8, height = 8)
plot(cv_fit$glmnet.fit,
     xvar = "lambda",
     main = "",
     col = colorRampPalette(c("steelblue", "firebrick"))(20),
     lwd = 1.5,
     xlab = expression(log(lambda)),
     ylab = "Coefficient Values")
abline(v = log(cv_fit$lambda.min), col = palette[1], lty = 2, lwd = 2)
abline(v = log(cv_fit$lambda.1se), col = palette[3], lty = 3, lwd = 2)
legend("topright",
       legend = c("Optimal lambda", "1-SE lambda"),
       col = palette[c(1,3)], lty = c(2,3), lwd = 2,
       cex = 0.9, bty = "n")
dev.off()
tiff("Coefficient_Path_Optimized.tif", width = 8, height = 8, units = "in", res = 300, compression = "lzw")
plot(cv_fit$glmnet.fit,
     xvar = "lambda",
     main = "",
     col = colorRampPalette(c("steelblue", "firebrick"))(20),
     lwd = 1.5,
     xlab = expression(log(lambda)),
     ylab = "Coefficient Values")
abline(v = log(cv_fit$lambda.min), col = palette[1], lty = 2, lwd = 2)
abline(v = log(cv_fit$lambda.1se), col = palette[3], lty = 3, lwd = 2)
legend("topright",
       legend = c("Optimal lambda", "1-SE lambda"),
       col = palette[c(1,3)], lty = c(2,3), lwd = 2,
       cex = 0.9, bty = "n")
dev.off()
roc_obj <- roc(y, risk_scores_vector)
pdf("ROC_Curve_Optimized.pdf", width = 8, height = 8)
plot(roc_obj,
     main = "ROC Curve for LASSO Model",
     col = palette[4],
     lwd = 3,
     print.auc = TRUE,
     print.auc.x = 0.6,
     print.auc.y = 0.2,
     print.auc.cex = 1.2,
     grid = TRUE,
     grid.col = "gray90",
     legacy.axes = TRUE)
legend("bottomright",
       legend = paste0("AUC = ", round(auc(roc_obj), 3)),
       col = palette[4], lwd = 3,
       cex = 1.1, bty = "n")
dev.off()
tiff("ROC_Curve_Optimized.tif", width = 8, height = 8, units = "in", res = 300, compression = "lzw")
plot(roc_obj,
     main = "ROC Curve for LASSO Model",
     col = palette[4],
     lwd = 3,
     print.auc = TRUE,
     print.auc.x = 0.6,
     print.auc.y = 0.2,
     print.auc.cex = 1.2,
     grid = TRUE,
     grid.col = "gray90",
     legacy.axes = TRUE)
legend("bottomright",
       legend = paste0("AUC = ", round(auc(roc_obj), 3)),
       col = palette[4], lwd = 3,
       cex = 1.1, bty = "n")
dev.off()
gene_order <- order(abs(selected_genes), decreasing = TRUE)
pdf("Feature_Importance.pdf", width = 8, height = 8)
barplot(selected_genes[gene_order],
        horiz = TRUE,
        las = 1,
        col = ifelse(selected_genes[gene_order] > 0, palette[2], palette[4]),
        border = NA,
        main = "Selected Genes by LASSO Model",
        xlab = "Coefficient Value",
        cex.names = 0.7,
        space = 0.5)
abline(v = 0, col = "gray30")
legend("topright",
       legend = c("Positive Association", "Negative Association"),
       fill = palette[c(2,4)],
       cex = 0.9, bty = "n")
dev.off()
tiff("Feature_Importance.tif", width = 8, height = 8, units = "in", res = 300, compression = "lzw")
barplot(selected_genes[gene_order],
        horiz = TRUE,
        las = 1,
        col = ifelse(selected_genes[gene_order] > 0, palette[2], palette[4]),
        border = NA,
        main = "Selected Genes by LASSO Model",
        xlab = "Coefficient Value",
        cex.names = 0.7,
        space = 0.5)
abline(v = 0, col = "gray30")
legend("topright",
       legend = c("Positive Association", "Negative Association"),
       fill = palette[c(2,4)],
       cex = 0.9, bty = "n")
dev.off()
selected_genes_df <- data.frame(
  Gene = names(selected_genes),
  Coefficient = selected_genes,
  row.names = NULL
)
write.csv(selected_genes_df[order(-abs(selected_genes_df$Coefficient)), ],
          "LASSO_Selected_Genes.csv", row.names = FALSE)
cat("=== Model Summary ===\n")
cat("- Optimal lambda (max AUC):", cv_fit$lambda.min, "\n")
cat("- 1-SE rule lambda:", cv_fit$lambda.1se, "\n")
cat("- Number of filtered genes:", length(selected_genes), "\n")
cat("- Model AUC:", round(auc(roc_obj), 3), "\n")
cat("- Risk score and filtered genes saved in: GSE25066-LASSO2.csv\n\n")
saveRDS(cv_fit, file = "lasso_model.rds")
