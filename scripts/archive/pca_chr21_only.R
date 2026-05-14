# pca_chr21_only.R
# PCA using only chr21 genes to see if karyotype separates

library(tidyverse)
library(DESeq2)

# Load data
count_data <- read_csv("data/processed/count_matrix.csv", show_col_types = FALSE)
metadata <- read_csv("data/processed/sample_metadata.csv", show_col_types = FALSE)
gene_annotations <- read_csv("data/processed/gene_annotations.csv", show_col_types = FALSE)

# Get count matrix
count_matrix <- count_data %>%
  select(-Gene_name, -Chr) %>%
  column_to_rownames("EnsemblID") %>%
  as.matrix()

# Round to integers
count_matrix <- round(count_matrix)

# Get chr21 genes only
chr21_genes <- gene_annotations %>%
  filter(Chr == "chr21") %>%
  pull(EnsemblID)

chr21_count_matrix <- count_matrix[chr21_genes, , drop = FALSE]

cat(sprintf("Chr21 genes: %d\n", nrow(chr21_count_matrix)))
cat(sprintf("Samples: %d\n", ncol(chr21_count_matrix)))

# Prepare colData
metadata <- metadata %>%
  filter(LabID %in% colnames(chr21_count_matrix)) %>%
  arrange(match(LabID, colnames(chr21_count_matrix)))

col_data <- DataFrame(
  sample_id = metadata$LabID,
  karyotype = factor(metadata$Karyotype, levels = c("Control", "T21")),
  sex = factor(metadata$Sex),
  age = metadata$Age_at_visit
)
rownames(col_data) <- metadata$LabID

# Create DESeqDataSet with chr21 genes only
dds_chr21 <- DESeqDataSetFromMatrix(
  countData = chr21_count_matrix,
  colData = col_data,
  design = ~ karyotype
)

# Filter lowly expressed genes
dds_chr21 <- dds_chr21[rowSums(counts(dds_chr21)) > 10, ]
cat(sprintf("Chr21 genes after filtering: %d\n\n", nrow(dds_chr21)))

# Variance stabilizing transformation
# Use varianceStabilizingTransformation for small gene sets
vsd <- varianceStabilizingTransformation(dds_chr21, blind = FALSE)

# Get PCA data
pca_data <- plotPCA(vsd, intgroup = "karyotype", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Create PCA plot
pdf("results/figures/pca_chr21_only.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

# Plot 1: All genes PCA (for comparison)
dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ karyotype
)
dds_all <- dds_all[rowSums(counts(dds_all)) > 10, ]
vsd_all <- vst(dds_all, blind = FALSE)
pca_all <- plotPCA(vsd_all, intgroup = "karyotype", returnData = TRUE)
percent_var_all <- round(100 * attr(pca_all, "percentVar"))

plot(pca_all$PC1, pca_all$PC2,
     col = ifelse(pca_all$karyotype == "T21", "red", "blue"),
     pch = 19, cex = 1.5, main = "PCA: All Genes",
     xlab = paste0("PC1: ", percent_var_all[1], "% variance"),
     ylab = paste0("PC2: ", percent_var_all[2], "% variance"))
legend("topright", legend = c("Control", "T21"),
       col = c("blue", "red"), pch = 19, cex = 1.2)

# Plot 2: Chr21 genes only
plot(pca_data$PC1, pca_data$PC2,
     col = ifelse(pca_data$karyotype == "T21", "red", "blue"),
     pch = 19, cex = 1.5, main = "PCA: Chr21 Genes Only",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"))
legend("topright", legend = c("Control", "T21"),
       col = c("blue", "red"), pch = 19, cex = 1.2)

dev.off()

cat("Saved: results/figures/pca_chr21_only.pdf\n\n")

# Print statistics
cat("=== PCA Statistics ===\n\n")

cat("All genes:\n")
cat(sprintf("  PC1 variance: %.1f%%\n", percent_var_all[1]))
cat(sprintf("  PC2 variance: %.1f%%\n", percent_var_all[2]))
cat(sprintf("  PC1 range: [%.2f, %.2f]\n", min(pca_all$PC1), max(pca_all$PC1)))
cat(sprintf("  PC2 range: [%.2f, %.2f]\n\n", min(pca_all$PC2), max(pca_all$PC2)))

cat("Chr21 genes only:\n")
cat(sprintf("  PC1 variance: %.1f%%\n", percent_var[1]))
cat(sprintf("  PC2 variance: %.1f%%\n", percent_var[2]))
cat(sprintf("  PC1 range: [%.2f, %.2f]\n", min(pca_data$PC1), max(pca_data$PC1)))
cat(sprintf("  PC2 range: [%.2f, %.2f]\n\n", min(pca_data$PC2), max(pca_data$PC2)))

# Check if PC1 separates by karyotype
cat("PC1 separation by karyotype:\n")
pc1_t21 <- pca_data$PC1[pca_data$karyotype == "T21"]
pc1_control <- pca_data$PC1[pca_data$karyotype == "Control"]

cat(sprintf("  Control PC1: mean=%.2f, sd=%.2f\n",
            mean(pc1_control), sd(pc1_control)))
cat(sprintf("  T21 PC1: mean=%.2f, sd=%.2f\n",
            mean(pc1_t21), sd(pc1_t21)))
cat(sprintf("  Difference: %.2f\n", mean(pc1_t21) - mean(pc1_control)))

# T-test
t_result <- t.test(pc1_t21, pc1_control)
cat(sprintf("  T-test p-value: %.2e\n", t_result$p.value))

if (t_result$p.value < 0.05) {
  cat("  *** PC1 SIGNIFICANTLY separates by karyotype! ***\n")
} else {
  cat("  PC1 does NOT significantly separate by karyotype\n")
}
