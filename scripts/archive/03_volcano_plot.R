# 03_volcano_plot.R
#
# Purpose: Create volcano plot of differential expression results
#          Highlight chromosome 21 genes and significant DE genes
#
# Inputs:
#   - results/tables/deseq2_all_genes_no_ploidy_norm.csv (raw FC for volcano)
#   - results/tables/chr21_genes_categorized.csv
#
# Outputs:
#   - results/figures/volcano_plot.pdf
#   - results/figures/volcano_plot_chr21_focus.pdf
#
# Date: 2025-11-11

# Load required libraries
library(tidyverse)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

cat("=== T21-eQTL Analysis: Volcano Plot ===\n\n")

# =============================================================================
# STEP 1: Load results
# =============================================================================

cat("Step 1: Loading DESeq2 results...\n")

if (!file.exists("results/tables/deseq2_all_genes_no_ploidy_norm.csv")) {
  stop("DESeq2 results not found. Run 01_deseq2_analysis.R first.")
}

# Use raw results (no ploidy normalization) for volcano plot
# This shows the actual fold changes observed in the data
results_all <- read_csv("results/tables/deseq2_all_genes_no_ploidy_norm.csv",
                        show_col_types = FALSE)

chr21_categorized <- read_csv("results/tables/chr21_genes_categorized.csv",
                              show_col_types = FALSE)

cat(sprintf("  Loaded %d genes total\n", nrow(results_all)))
cat(sprintf("  Chr21 genes: %d\n", nrow(chr21_categorized)))

# =============================================================================
# STEP 2: Prepare data for plotting
# =============================================================================

cat("\nStep 2: Preparing data for plotting...\n")

# Add significance and chromosome labels
plot_data <- results_all %>%
  mutate(
    # Negative log10 p-value for y-axis
    neg_log10_padj = -log10(padj),
    # Significance categories
    significance = case_when(
      is.na(padj) ~ "Not tested",
      padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant (FC>=2)",
      padj < 0.05 ~ "Significant (FC<2)",
      TRUE ~ "Not significant"
    ),
    # Chromosome categories
    chr_category = case_when(
      Chr == "chr21" ~ "Chr21",
      !is.na(Chr) ~ "Other chromosomes",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Remove genes with NA padj for cleaner plot
  filter(!is.na(padj))

cat(sprintf("  Genes after filtering: %d\n", nrow(plot_data)))

# Summary of significance
sig_summary <- plot_data %>%
  count(significance)
cat("\n  Significance breakdown:\n")
print(as.data.frame(sig_summary))

# =============================================================================
# STEP 3: Identify genes to label
# =============================================================================

cat("\nStep 3: Identifying genes to label...\n")

# Label criteria:
# 1. Top significant chr21 genes
# 2. Top significant non-chr21 genes
# 3. Chr21 genes with interesting patterns

# Top chr21 genes
top_chr21 <- plot_data %>%
  filter(Chr == "chr21", padj < 0.05) %>%
  arrange(padj) %>%
  head(10)

# Top other genes
top_other <- plot_data %>%
  filter(Chr != "chr21", padj < 0.05) %>%
  arrange(padj) %>%
  head(5)

# Combine genes to label
genes_to_label <- bind_rows(top_chr21, top_other) %>%
  distinct(EnsemblID, .keep_all = TRUE)

cat(sprintf("  Genes to label: %d\n", nrow(genes_to_label)))

# =============================================================================
# STEP 4: Create volcano plot - all chromosomes
# =============================================================================

cat("\nStep 4: Creating volcano plot (all chromosomes)...\n")

# Color palette
colors <- c(
  "Chr21" = "#E41A1C",
  "Other chromosomes" = "gray60"
)

# Create plot
p1 <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
  # Background points (other chromosomes)
  geom_point(data = filter(plot_data, chr_category == "Other chromosomes"),
             aes(color = chr_category),
             alpha = 0.3, size = 1) +
  # Chr21 points on top
  geom_point(data = filter(plot_data, chr_category == "Chr21"),
             aes(color = chr_category),
             alpha = 0.7, size = 2) +
  # Reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "black", alpha = 0.3) +
  # Labels for significant genes
  geom_text_repel(data = genes_to_label,
                  aes(label = Gene_name),
                  size = 3, max.overlaps = 20,
                  box.padding = 0.5) +
  # Scales and labels
  scale_color_manual(values = colors) +
  labs(
    title = "Differential Expression: T21 vs Control",
    subtitle = "Trisomy-aware DESeq2 with ploidy normalization",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)",
    color = "Chromosome"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )

# Save plot
ggsave("results/figures/volcano_plot.pdf", p1,
       width = 10, height = 8)
cat("  Saved: results/figures/volcano_plot.pdf\n")

# =============================================================================
# STEP 5: Create volcano plot - chr21 focus
# =============================================================================

cat("\nStep 5: Creating volcano plot (chr21 focus)...\n")

# Prepare chr21 data with categories
chr21_plot_data <- plot_data %>%
  filter(Chr == "chr21") %>%
  left_join(chr21_categorized %>%
              select(EnsemblID, category),
            by = "EnsemblID")

# Color by category
category_colors <- c(
  "Expected Dosage (>=1.5 FC)" = "#377EB8",
  "High Genomic Repeats" = "#E41A1C",
  "Low Gene Expression" = "#984EA3",
  "Not Differentially Expressed" = "#999999",
  "Differentially Expressed (<1.5 FC)" = "#FF7F00"
)

p2 <- ggplot(chr21_plot_data,
             aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = category), alpha = 0.7, size = 3) +
  # Reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "red", alpha = 0.5, linewidth = 1) +
  # Label all significant chr21 genes
  geom_text_repel(data = filter(chr21_plot_data, padj < 0.05),
                  aes(label = Gene_name),
                  size = 3, max.overlaps = 30) +
  # Scales and labels
  scale_color_manual(values = category_colors,
                     name = "Category") +
  labs(
    title = "Chromosome 21 Genes: T21 vs Control",
    subtitle = paste("After ploidy normalization, expected log2FC = 0",
                     "(genes at log2FC >= 0 show expected 1.5-fold change)"),
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9),
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

# Save plot
ggsave("results/figures/volcano_plot_chr21_focus.pdf", p2,
       width = 12, height = 8)
cat("  Saved: results/figures/volcano_plot_chr21_focus.pdf\n")

# =============================================================================
# STEP 6: Summary statistics
# =============================================================================

cat("\n=== Volcano Plot Summary ===\n")

# Overall statistics
cat(sprintf("Total genes plotted: %d\n", nrow(plot_data)))
cat(sprintf("Chr21 genes: %d\n", sum(plot_data$Chr == "chr21")))
cat(sprintf("Significant genes (padj < 0.05): %d\n",
            sum(plot_data$padj < 0.05)))

# Chr21 specific
chr21_sig <- sum(plot_data$Chr == "chr21" & plot_data$padj < 0.05)
chr21_total <- sum(plot_data$Chr == "chr21")
cat(sprintf("\nChr21 significant: %d / %d (%.1f%%)\n",
            chr21_sig, chr21_total, 100 * chr21_sig / chr21_total))

# Fold change distribution for chr21
chr21_fc_pos <- sum(plot_data$Chr == "chr21" &
                    plot_data$log2FoldChange >= 0, na.rm = TRUE)
chr21_fc_neg <- sum(plot_data$Chr == "chr21" &
                    plot_data$log2FoldChange < 0, na.rm = TRUE)
cat(sprintf("Chr21 genes at expected dosage (log2FC >= 0): %d (%.1f%%)\n",
            chr21_fc_pos, 100 * chr21_fc_pos / chr21_total))
cat(sprintf("Chr21 genes below expected (log2FC < 0): %d (%.1f%%)\n",
            chr21_fc_neg, 100 * chr21_fc_neg / chr21_total))

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/figures/volcano_plot_session_info.txt")

cat("\n=== Volcano Plot Complete ===\n")
cat("Next step: Run 04_alluvial_plot.R\n\n")
