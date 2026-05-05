# diagnostic_check.R
# Quick diagnostic to check if chr21 genes are actually elevated in T21

library(tidyverse)

# Load data
count_data <- read_csv("data/processed/count_matrix.csv", show_col_types = FALSE)
metadata <- read_csv("data/processed/sample_metadata.csv", show_col_types = FALSE)
gene_annotations <- read_csv("data/processed/gene_annotations.csv", show_col_types = FALSE)

# Get count matrix
count_matrix <- count_data %>%
  select(-Gene_name, -Chr) %>%
  column_to_rownames("EnsemblID") %>%
  as.matrix()

# Get chr21 genes
chr21_genes <- gene_annotations %>%
  filter(Chr == "chr21") %>%
  pull(EnsemblID)

# Get T21 and Control samples
t21_samples <- metadata %>% filter(Karyotype == "T21") %>% pull(LabID)
control_samples <- metadata %>% filter(Karyotype == "Control") %>% pull(LabID)

cat(sprintf("Samples: %d T21, %d Control\n", length(t21_samples), length(control_samples)))
cat(sprintf("Chr21 genes: %d\n\n", length(chr21_genes)))

# Calculate mean counts for chr21 genes
chr21_counts <- count_matrix[chr21_genes, , drop = FALSE]

# Remove lowly expressed genes
chr21_expressed <- chr21_counts[rowMeans(chr21_counts) > 10, , drop = FALSE]
cat(sprintf("Chr21 genes with mean > 10: %d\n\n", nrow(chr21_expressed)))

# Calculate mean expression in each group
t21_mean <- rowMeans(chr21_expressed[, t21_samples, drop = FALSE])
control_mean <- rowMeans(chr21_expressed[, control_samples, drop = FALSE])

# Calculate raw fold change
raw_fc <- t21_mean / control_mean
log2_fc <- log2(raw_fc)

cat("Chr21 gene raw fold changes (T21 / Control):\n")
cat(sprintf("  Median FC: %.3f\n", median(raw_fc, na.rm = TRUE)))
cat(sprintf("  Median log2FC: %.3f (expected: %.3f for 1.5x)\n",
            median(log2_fc, na.rm = TRUE), log2(1.5)))
cat(sprintf("  Mean FC: %.3f\n", mean(raw_fc, na.rm = TRUE)))
cat(sprintf("  Mean log2FC: %.3f\n\n", mean(log2_fc, na.rm = TRUE)))

cat("Distribution of fold changes:\n")
print(summary(raw_fc))

cat("\nNumber of genes by FC category:\n")
cat(sprintf("  FC >= 1.5: %d (%.1f%%)\n",
            sum(raw_fc >= 1.5, na.rm = TRUE),
            100 * sum(raw_fc >= 1.5, na.rm = TRUE) / length(raw_fc)))
cat(sprintf("  FC < 1.5: %d (%.1f%%)\n",
            sum(raw_fc < 1.5, na.rm = TRUE),
            100 * sum(raw_fc < 1.5, na.rm = TRUE) / length(raw_fc)))
cat(sprintf("  FC < 1.0: %d (%.1f%%)\n",
            sum(raw_fc < 1.0, na.rm = TRUE),
            100 * sum(raw_fc < 1.0, na.rm = TRUE) / length(raw_fc)))

# Check a few example genes
cat("\nExample chr21 genes:\n")
example_genes <- head(chr21_expressed, 10)
for (i in 1:min(10, nrow(example_genes))) {
  gene_name <- gene_annotations$Gene_name[gene_annotations$EnsemblID == rownames(example_genes)[i]]
  cat(sprintf("  %s: T21=%.1f, Control=%.1f, FC=%.2f\n",
              gene_name,
              t21_mean[i],
              control_mean[i],
              raw_fc[i]))
}

# Check overall library sizes
cat("\n\nLibrary size check:\n")
total_counts_t21 <- colSums(count_matrix[, t21_samples])
total_counts_control <- colSums(count_matrix[, control_samples])
cat(sprintf("  T21 median total counts: %.0f\n", median(total_counts_t21)))
cat(sprintf("  Control median total counts: %.0f\n", median(total_counts_control)))
cat(sprintf("  Ratio: %.3f\n", median(total_counts_t21) / median(total_counts_control)))
