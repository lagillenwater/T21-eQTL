# 08_expected_dosage_eqtl.R
#
# Purpose: Analyze chr21 genes with >= 1.5 FC that remain significant after
#          ploidy correction. These genes deviate from expected dosage even
#          after accounting for trisomy. Cross-reference with eQTL data.
#
# Research question: Among genes showing expected 1.5x dosage increase,
#                    are any still differentially expressed after ploidy
#                    normalization? If so, can eQTLs explain this deviation?
#
# Inputs:
#   - results/tables/chr21_genes_categorized.csv
#   - data/Whole_Blood.v10.eQTLs.signif_pairs.parquet (GTEx v10)
#
# Outputs:
#   - results/tables/expected_dosage_significant.csv
#   - results/tables/expected_dosage_eqtl_matches.csv
#   - results/figures/expected_dosage_eqtl_summary.pdf
#
# Author: Claude Code
# Date: 2025-12-17

# Load required libraries
library(tidyverse)
library(arrow)

# Set seed for reproducibility
set.seed(42)

cat("=== Expected Dosage Gene Analysis: eQTL Cross-Reference ===\n\n")

# =============================================================================
# STEP 1: Load categorized chr21 genes
# =============================================================================

cat("Step 1: Loading categorized chr21 genes...\n")

if (!file.exists("results/tables/chr21_genes_categorized.csv")) {
  stop("Categorized data not found. Run 02_categorize_genes.R first.")
}

chr21_categorized <- read_csv("results/tables/chr21_genes_categorized.csv",
                              show_col_types = FALSE)

# =============================================================================
# STEP 2: Identify genes with >= 1.5 FC that are still significant after
#         ploidy correction
# =============================================================================

cat("\nStep 2: Identifying expected dosage genes significant after ploidy correction...\n")

# Get genes with >= 1.5 FC (expected dosage)
expected_dosage_genes <- chr21_categorized %>%
  filter(category == "Expected Dosage (>=1.5 FC)")

n_expected_total <- nrow(expected_dosage_genes)
cat(sprintf("  Total expected dosage genes (>=1.5 FC): %d\n", n_expected_total))

# Filter for those still significant after ploidy normalization
# These genes deviate from expected even after accounting for trisomy
expected_dosage_significant <- expected_dosage_genes %>%
  filter(norm_padj < 0.05) %>%
  arrange(norm_padj) %>%
  mutate(
    EnsemblID_clean = str_remove(EnsemblID, "\\..*"),
    # Classify direction of deviation
    deviation_direction = case_when(
      norm_log2FC > 0.1 ~ "Higher than expected",
      norm_log2FC < -0.1 ~ "Lower than expected",
      TRUE ~ "Near expected"
    )
  )

n_significant <- nrow(expected_dosage_significant)
cat(sprintf("  Significant after ploidy correction (padj < 0.05): %d\n",
            n_significant))
cat(sprintf("  Percentage: %.1f%%\n\n", 100 * n_significant / n_expected_total))

if (n_significant == 0) {
  cat("=== No expected dosage genes remain significant after ploidy correction ===\n")
  cat("This indicates that ploidy normalization fully accounts for their expression.\n")

  # Create empty output
  empty_results <- data.frame(
    Gene_name = character(),
    EnsemblID = character(),
    category = character()
  )
  write_csv(empty_results, "results/tables/expected_dosage_significant.csv")
  quit(save = "no")
}

# Summarize deviation directions
cat("  Deviation from expected dosage:\n")
deviation_summary <- expected_dosage_significant %>%
  count(deviation_direction) %>%
  mutate(percentage = round(100 * n / sum(n), 1))
print(as.data.frame(deviation_summary))

cat("\n  Top genes by significance:\n")
top_genes <- expected_dosage_significant %>%
  select(Gene_name, raw_log2FC, norm_log2FC, raw_padj, norm_padj,
         baseMean, deviation_direction) %>%
  head(15)
print(as.data.frame(top_genes))

# Save intermediate results
write_csv(expected_dosage_significant,
          "results/tables/expected_dosage_significant.csv")
cat("\n  Saved: results/tables/expected_dosage_significant.csv\n")

# =============================================================================
# STEP 3: Load GTEx eQTL data
# =============================================================================

cat("\nStep 3: Loading GTEx eQTL data...\n")

gtex_file <- "data/Whole_Blood.v10.eQTLs.signif_pairs.parquet"

if (!file.exists(gtex_file)) {
  cat("\n*** GTEx eQTL data not found ***\n")
  cat("Please download from GTEx Portal:\n")
  cat("  https://gtexportal.org/home/datasets\n")
  cat("  File: Whole_Blood.v10.eQTLs.signif_pairs.parquet\n")
  cat("  Save to: data/ directory\n\n")

  # Create template for manual lookup
  eqtl_template <- expected_dosage_significant %>%
    select(Gene_name, EnsemblID, raw_log2FC, norm_log2FC, baseMean,
           norm_padj, deviation_direction) %>%
    mutate(
      eQTL_found = NA,
      eQTL_count = NA,
      notes = "Manual lookup required - GTEx file not found"
    )

  write_csv(eqtl_template,
            "results/tables/expected_dosage_eqtl_template.csv")
  cat("Created template: results/tables/expected_dosage_eqtl_template.csv\n")
  quit(save = "no")
}

cat("  Reading GTEx v10 parquet file...\n")
gtex_eqtls <- read_parquet(gtex_file)
cat(sprintf("  Loaded %d total eQTL associations\n", nrow(gtex_eqtls)))

# Check column names
cat("  GTEx columns: ", paste(colnames(gtex_eqtls)[1:10], collapse = ", "), "...\n")

# =============================================================================
# STEP 4: Filter GTEx data for our genes of interest
# =============================================================================

cat("\nStep 4: Filtering GTEx eQTLs for expected dosage genes...\n")

# Clean gene IDs for matching
gtex_filtered <- gtex_eqtls %>%
  mutate(gene_id_clean = str_remove(gene_id, "\\..*")) %>%
  filter(gene_id_clean %in% expected_dosage_significant$EnsemblID_clean)

n_eqtl_associations <- nrow(gtex_filtered)
n_genes_with_eqtls <- length(unique(gtex_filtered$gene_id_clean))

cat(sprintf("  eQTL associations found: %d\n", n_eqtl_associations))
cat(sprintf("  Genes with at least one eQTL: %d / %d\n",
            n_genes_with_eqtls, n_significant))

# =============================================================================
# STEP 5: Summarize eQTL data per gene
# =============================================================================

cat("\nStep 5: Summarizing eQTLs per gene...\n")

# Aggregate eQTL information per gene
eqtl_summary <- gtex_filtered %>%
  group_by(gene_id_clean) %>%
  summarize(
    eQTL_count = n(),
    n_positive_slope = sum(slope > 0, na.rm = TRUE),
    n_negative_slope = sum(slope < 0, na.rm = TRUE),
    mean_slope = mean(slope, na.rm = TRUE),
    min_pval = min(pval_nominal, na.rm = TRUE),
    top_variant = variant_id[which.min(pval_nominal)],
    .groups = "drop"
  ) %>%
  mutate(
    predominant_direction = case_when(
      n_positive_slope > n_negative_slope ~ "Positive (higher expression)",
      n_negative_slope > n_positive_slope ~ "Negative (lower expression)",
      TRUE ~ "Mixed"
    )
  )

# =============================================================================
# STEP 6: Match genes with eQTL results
# =============================================================================

cat("\nStep 6: Matching genes with eQTL data...\n")

eqtl_matches <- expected_dosage_significant %>%
  left_join(eqtl_summary, by = c("EnsemblID_clean" = "gene_id_clean")) %>%
  mutate(
    eQTL_found = !is.na(eQTL_count),
    # Interpret whether eQTL explains the deviation
    eqtl_interpretation = case_when(
      !eQTL_found ~ "No eQTL found",
      # Higher expression + positive eQTL = eQTL may explain
      deviation_direction == "Higher than expected" &
        predominant_direction == "Positive (higher expression)" ~
        "eQTL supports higher expression",
      # Lower expression + negative eQTL = eQTL may explain
      deviation_direction == "Lower than expected" &
        predominant_direction == "Negative (lower expression)" ~
        "eQTL supports lower expression",
      # Mismatch between deviation and eQTL direction
      TRUE ~ "eQTL direction does not match deviation"
    )
  )

# Summary statistics
cat("\n  eQTL interpretation summary:\n")
interpretation_summary <- eqtl_matches %>%
  count(eqtl_interpretation) %>%
  mutate(percentage = round(100 * n / sum(n), 1))
print(as.data.frame(interpretation_summary))

# =============================================================================
# STEP 7: Detailed results table
# =============================================================================

cat("\nStep 7: Creating detailed results...\n")

results_table <- eqtl_matches %>%
  select(
    Gene_name, EnsemblID, baseMean,
    raw_log2FC, norm_log2FC, norm_padj,
    deviation_direction,
    eQTL_found, eQTL_count, mean_slope, predominant_direction,
    top_variant, eqtl_interpretation
  ) %>%
  arrange(norm_padj)

write_csv(results_table, "results/tables/expected_dosage_eqtl_matches.csv")
cat("  Saved: results/tables/expected_dosage_eqtl_matches.csv\n")

# Print detailed results
cat("\n=== Detailed Gene Results ===\n\n")
for (i in seq_len(nrow(results_table))) {
  gene <- results_table[i, ]
  cat(sprintf("%d. %s\n", i, gene$Gene_name))
  cat(sprintf("   Raw log2FC: %.3f | Normalized log2FC: %.3f\n",
              gene$raw_log2FC, gene$norm_log2FC))
  cat(sprintf("   Ploidy-corrected padj: %.2e\n", gene$norm_padj))
  cat(sprintf("   Deviation: %s\n", gene$deviation_direction))
  if (gene$eQTL_found) {
    cat(sprintf("   eQTLs: %d associations, mean slope: %.3f (%s)\n",
                gene$eQTL_count, gene$mean_slope, gene$predominant_direction))
    cat(sprintf("   Top variant: %s\n", gene$top_variant))
  } else {
    cat("   eQTLs: None found in GTEx whole blood\n")
  }
  cat(sprintf("   Interpretation: %s\n\n", gene$eqtl_interpretation))
}

# =============================================================================
# STEP 8: Create visualizations
# =============================================================================

cat("Step 8: Creating visualizations...\n")

# Plot 1: Bar chart of eQTL interpretation
p_interpretation <- ggplot(
  interpretation_summary,
  aes(x = reorder(eqtl_interpretation, -n), y = n, fill = eqtl_interpretation)
) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n, percentage)),
            vjust = -0.5, size = 4) +
  scale_fill_manual(
    values = c(
      "eQTL supports higher expression" = "#4DAF4A",
      "eQTL supports lower expression" = "#377EB8",
      "eQTL direction does not match deviation" = "#FF7F00",
      "No eQTL found" = "#999999"
    )
  ) +
  labs(
    title = "eQTL Analysis: Expected Dosage Genes Significant After Ploidy Correction",
    subtitle = sprintf("n = %d genes with >=1.5 FC that remain significant (padj < 0.05)",
                       n_significant),
    x = "",
    y = "Number of Genes"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    legend.position = "none"
  ) +
  ylim(0, max(interpretation_summary$n) * 1.15)

# Plot 2: Scatter of norm_log2FC vs eQTL mean slope
genes_with_eqtls <- eqtl_matches %>% filter(eQTL_found)

if (nrow(genes_with_eqtls) > 0) {
  p_scatter <- ggplot(
    genes_with_eqtls,
    aes(x = mean_slope, y = norm_log2FC, size = eQTL_count)
  ) +
    geom_point(aes(color = eqtl_interpretation), alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = Gene_name), vjust = -1, hjust = 0.5, size = 3) +
    scale_color_manual(
      values = c(
        "eQTL supports higher expression" = "#4DAF4A",
        "eQTL supports lower expression" = "#377EB8",
        "eQTL direction does not match deviation" = "#FF7F00"
      )
    ) +
    labs(
      title = "Expression Deviation vs eQTL Effect Direction",
      subtitle = "Expected dosage genes with known eQTLs",
      x = "Mean eQTL Slope (GTEx Whole Blood)",
      y = "Normalized log2 Fold Change\n(deviation from expected)",
      color = "Interpretation",
      size = "# eQTLs"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "right"
    )
} else {
  p_scatter <- NULL
}

# Plot 3: Distribution of normalized log2FC by eQTL status
p_boxplot <- ggplot(
  eqtl_matches,
  aes(x = eQTL_found, y = norm_log2FC, fill = eQTL_found)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("FALSE" = "No eQTL", "TRUE" = "Has eQTL")) +
  scale_fill_manual(values = c("FALSE" = "#999999", "TRUE" = "#377EB8")) +
  labs(
    title = "Expression Deviation by eQTL Status",
    subtitle = "Expected dosage genes significant after ploidy correction",
    x = "eQTL Status (GTEx Whole Blood)",
    y = "Normalized log2 Fold Change\n(deviation from expected)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

# Save plots
pdf("results/figures/expected_dosage_eqtl_summary.pdf", width = 12, height = 10)
print(p_interpretation)
if (!is.null(p_scatter)) print(p_scatter)
print(p_boxplot)
dev.off()
cat("  Saved: results/figures/expected_dosage_eqtl_summary.pdf\n")

# =============================================================================
# STEP 9: Final summary
# =============================================================================

cat("\n=== FINAL SUMMARY ===\n\n")

cat(sprintf("Total chr21 genes analyzed: %d\n", nrow(chr21_categorized)))
cat(sprintf("Expected dosage genes (>=1.5 FC): %d\n", n_expected_total))
cat(sprintf("Significant after ploidy correction: %d (%.1f%%)\n\n",
            n_significant, 100 * n_significant / n_expected_total))

# Genes with and without eQTLs
n_with_eqtl <- sum(eqtl_matches$eQTL_found)
n_without_eqtl <- sum(!eqtl_matches$eQTL_found)

cat("eQTL status:\n")
cat(sprintf("  With known eQTLs: %d (%.1f%%)\n",
            n_with_eqtl, 100 * n_with_eqtl / n_significant))
cat(sprintf("  No known eQTLs: %d (%.1f%%)\n\n",
            n_without_eqtl, 100 * n_without_eqtl / n_significant))

# Interpretation breakdown
cat("Interpretation:\n")
for (i in seq_len(nrow(interpretation_summary))) {
  row <- interpretation_summary[i, ]
  cat(sprintf("  %s: %d (%.1f%%)\n",
              row$eqtl_interpretation, row$n, row$percentage))
}

# Biological interpretation
cat("\n=== BIOLOGICAL INTERPRETATION ===\n\n")

n_eqtl_supported <- sum(grepl("supports", eqtl_matches$eqtl_interpretation))

cat("These genes show > 1.5-fold change (expected for chr21 in T21) but\n")
cat("remain significantly different from expected even after ploidy correction.\n\n")

if (n_eqtl_supported > 0) {
  cat(sprintf("%d/%d genes have eQTL effects that align with their deviation\n",
              n_eqtl_supported, n_significant))
  cat("direction, suggesting genetic variants contribute to expression\n")
  cat("variation beyond copy number.\n\n")

  supported_genes <- eqtl_matches %>%
    filter(grepl("supports", eqtl_interpretation)) %>%
    pull(Gene_name)
  cat("Genes with supporting eQTLs: ", paste(supported_genes, collapse = ", "), "\n\n")
}

unsupported_genes <- eqtl_matches %>%
  filter(eqtl_interpretation == "No eQTL found" |
         eqtl_interpretation == "eQTL direction does not match deviation") %>%
  pull(Gene_name)

if (length(unsupported_genes) > 0) {
  cat(sprintf("%d genes lack eQTL explanation for their deviation:\n",
              length(unsupported_genes)))
  cat("  ", paste(unsupported_genes, collapse = ", "), "\n\n")
  cat("These genes may be regulated by:\n")
  cat("  - Rare variants not captured by common eQTL studies\n")
  cat("  - Tissue-specific or context-dependent regulation\n")
  cat("  - Epigenetic mechanisms\n")
  cat("  - Trans-acting factors affected by trisomy\n")
}

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/tables/expected_dosage_eqtl_session_info.txt")

cat("\n=== Analysis Complete ===\n")
cat("All results saved to results/ directory.\n\n")
