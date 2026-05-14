# 05_eqtl_analysis.R
#
# Purpose: Cross-reference differentially expressed chr21 genes with eQTL databases
#          Determine if low expression can be explained by known genetic variants
#
# Inputs:
#   - results/tables/chr21_genes_categorized.csv
#   - data/gtex_whole_blood_eqtls.csv (user must download from GTEx Portal)
#
# Outputs:
#   - results/tables/eqtl_matches.csv
#   - results/tables/final_gene_classification.csv
#   - results/figures/eqtl_summary.pdf
#
# Date: 2025-11-11
#
# NOTE: This script requires GTEx eQTL data to be downloaded separately
#       from https://gtexportal.org/home/datasets

# Load required libraries
library(tidyverse)
library(arrow)  # For reading parquet files

# Set seed for reproducibility
set.seed(42)

cat("=== T21-eQTL Analysis: eQTL Cross-Reference ===\n\n")

# =============================================================================
# STEP 1: Load categorized chr21 genes
# =============================================================================

cat("Step 1: Loading categorized chr21 genes...\n")

if (!file.exists("results/tables/chr21_genes_categorized.csv")) {
  stop("Categorized data not found. Run 02_categorize_genes.R first.")
}

chr21_categorized <- read_csv("results/tables/chr21_genes_categorized.csv",
                              show_col_types = FALSE)

# Get differentially expressed genes with low fold change
# These are genes with raw FC < 1.5 that remain significant after ploidy normalization
de_genes_low_fc <- chr21_categorized %>%
  filter(category == "Expression/Ploidy Differentially Expressed") %>%
  arrange(norm_padj) %>%
  mutate(EnsemblID_clean = str_remove(EnsemblID, "\\..*"))

n_de_genes <- nrow(de_genes_low_fc)
cat(sprintf("  Total chr21 genes: %d\n", nrow(chr21_categorized)))
cat(sprintf("  DE genes with FC < 1.5: %d\n", n_de_genes))

if (n_de_genes == 0) {
  cat("\n=== No DE genes found for eQTL analysis ===\n")
  cat("This is consistent with NO dosage compensation in chr21.\n")
  cat("Analysis complete - no further eQTL lookup needed.\n\n")

  # Create empty output files
  empty_results <- data.frame(
    Gene_name = character(),
    EnsemblID = character(),
    eQTL_found = logical(),
    eQTL_count = integer(),
    eQTL_explains_expression = character(),
    notes = character()
  )

  write_csv(empty_results, "results/tables/eqtl_matches.csv")
  write_csv(chr21_categorized, "results/tables/final_gene_classification.csv")

  quit(save = "no")
}

cat("\n  Genes requiring eQTL analysis:\n")
de_summary <- de_genes_low_fc %>%
  select(Gene_name, raw_log2FC, norm_log2FC, baseMean, raw_padj, norm_padj)
print(as.data.frame(de_summary))

# =============================================================================
# STEP 2: Download and load GTEx eQTL data
# =============================================================================

cat("\nStep 2: Checking for GTEx eQTL data...\n")

# GTEx Whole Blood eQTL file (check for v10 parquet or v8 text file)
gtex_file_v10 <- "data/Whole_Blood.v10.eQTLs.signif_pairs.parquet"
gtex_file_v8 <- "data/gtex_whole_blood_eqtls.txt.gz"

# Check which file is available
if (file.exists(gtex_file_v10)) {
  cat("  Found GTEx v10 parquet file (local)\n")
  gtex_file <- gtex_file_v10
  gtex_format <- "parquet"
} else if (file.exists(gtex_file_v8)) {
  cat("  Found GTEx v8 text file (local)\n")
  gtex_file <- gtex_file_v8
  gtex_format <- "txt.gz"
} else {
  # Try to download v8 format
  cat("  GTEx eQTL data not found locally.\n")
  gtex_file <- gtex_file_v8
  gtex_format <- "txt.gz"

# Try multiple possible URLs (GTEx has moved files between versions)
gtex_urls <- c(
  # Current GTEx Portal (v8)
  "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",
  # Alternative GTEx storage
  "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",
  # GTEx download portal direct link
  "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/all_snp_gene_pairs/Whole_Blood.allpairs.txt.gz"
)

  # No file found - provide manual instructions
  cat("\n*** GTEx eQTL data not found ***\n")
  cat("\nPlease download manually:\n")
  cat("1. Visit: https://gtexportal.org/home/datasets\n")
  cat("2. Navigate to: Single-Tissue cis-eQTL Data\n")
  cat("3. Download one of:\n")
  cat("   - Whole_Blood.v10.eQTLs.signif_pairs.parquet (recommended, v10)\n")
  cat("   - Whole_Blood.v8.signif_variant_gene_pairs.txt.gz (v8)\n")
  cat("4. Save to: data/ directory\n")

  # Create template for manual lookup
  eqtl_template <- de_genes_low_fc %>%
    select(Gene_name, EnsemblID, raw_log2FC, norm_log2FC, baseMean,
           raw_padj, norm_padj) %>%
    mutate(
      eQTL_found = NA,
      eQTL_count = NA,
      eQTL_direction = NA,
      eQTL_explains_low_expression = NA,
      notes = "Manual lookup required - GTEx file not found"
    )

  write_csv(eqtl_template, "results/tables/eqtl_lookup_template.csv")
  cat("\nCreated template: results/tables/eqtl_lookup_template.csv\n")

  quit(save = "no")
}

# =============================================================================
# STEP 3: Load and process GTEx eQTL data
# =============================================================================

cat("\nStep 3: Loading GTEx eQTL data...\n")
cat("  This may take 1-2 minutes for large file...\n")

# Load GTEx data based on format
if (gtex_format == "parquet") {
  cat("  Reading parquet format...\n")
  gtex_eqtls <- read_parquet(gtex_file)
} else {
  cat("  Reading text format...\n")
  gtex_eqtls <- read_tsv(gtex_file, show_col_types = FALSE)
}

cat(sprintf("  Loaded %d total eQTL associations\n", nrow(gtex_eqtls)))

# Check required columns
required_cols <- c("gene_id", "variant_id", "pval_nominal", "slope")
if (!all(required_cols %in% colnames(gtex_eqtls))) {
  stop(paste("GTEx file missing required columns:",
             paste(setdiff(required_cols, colnames(gtex_eqtls)),
                   collapse = ", ")))
}

cat("  Filtering to chr21 genes and significant eQTLs...\n")

# Filter for chr21 genes only (starts with ENSG, on chr21)
# GTEx only includes significant variant-gene pairs (FDR < 0.05)
# We'll further filter to nominal p < 1e-5 for stronger associations
gtex_chr21 <- gtex_eqtls %>%
  # Remove Ensembl version numbers for matching
  mutate(gene_id_clean = str_remove(gene_id, "\\..*")) %>%
  # Keep only chr21 genes (filter to our DE gene list)
  filter(gene_id_clean %in% de_genes_low_fc$EnsemblID_clean) %>%
  # Further filter to strong associations
  filter(pval_nominal < 1e-5)

cat(sprintf("  Chr21 DE genes with significant eQTLs (p < 1e-5): %d associations\n",
            nrow(gtex_chr21)))

if (nrow(gtex_chr21) == 0) {
  cat("\n*** No significant eQTLs found for chr21 DE genes ***\n")
  cat("This suggests these genes do not have common eQTLs in whole blood.\n")

  # Create empty results
  gtex_significant <- gtex_chr21  # Empty data frame
} else {
  gtex_significant <- gtex_chr21

  # Show summary
  n_genes_with_eqtls <- length(unique(gtex_significant$gene_id_clean))
  cat(sprintf("  Number of chr21 DE genes with eQTLs: %d\n", n_genes_with_eqtls))
}

# =============================================================================
# STEP 4: Match DE genes with eQTLs
# =============================================================================

cat("\nStep 4: Matching DE genes with eQTLs...\n")

# Match with GTEx (EnsemblID_clean already created above)
eqtl_matches <- de_genes_low_fc %>%
  left_join(
    gtex_significant %>%
      group_by(gene_id_clean) %>%
      summarize(
        eQTL_count = n(),
        min_pval = min(pval_nominal),
        mean_slope = mean(slope),
        negative_slope_count = sum(slope < 0),
        .groups = "drop"
      ),
    by = c("EnsemblID_clean" = "gene_id_clean")
  ) %>%
  mutate(
    eQTL_found = !is.na(eQTL_count),
    # An eQTL explains low expression if:
    # 1. eQTL exists
    # 2. Majority of eQTLs have negative slope (lower expression)
    eQTL_explains_expression = case_when(
      !eQTL_found ~ "No eQTL found",
      negative_slope_count > (eQTL_count / 2) ~ "eQTL supports low expression",
      TRUE ~ "eQTL does not support low expression"
    )
  )

# Summary
n_with_eqtl <- sum(eqtl_matches$eQTL_found)
n_explained <- sum(eqtl_matches$eQTL_explains_expression ==
                   "eQTL supports low expression")
n_unexplained <- n_de_genes - n_explained

cat(sprintf("  Genes with eQTLs: %d / %d\n", n_with_eqtl, n_de_genes))
cat(sprintf("  Expression explained by eQTLs: %d\n", n_explained))
cat(sprintf("  Expression NOT explained by eQTLs: %d\n", n_unexplained))

# =============================================================================
# STEP 5: Create final classification
# =============================================================================

cat("\nStep 5: Creating final classification...\n")

# Add eQTL info to full categorization
final_classification <- chr21_categorized %>%
  left_join(
    eqtl_matches %>%
      select(EnsemblID, eQTL_found, eQTL_count, eQTL_explains_expression),
    by = "EnsemblID"
  ) %>%
  mutate(
    # Final category including eQTL info
    final_category = case_when(
      category != "Differentially Expressed (<1.5 FC)" ~ category,
      eQTL_explains_expression == "eQTL supports low expression" ~
        "DE: eQTL-supported",
      eQTL_explains_expression == "No eQTL found" ~
        "DE: No eQTL explanation",
      TRUE ~ "DE: eQTL does not support"
    )
  )

# =============================================================================
# STEP 6: Save results
# =============================================================================

cat("\nStep 6: Saving results...\n")

# Save eQTL matches
eqtl_results <- eqtl_matches %>%
  select(Gene_name, EnsemblID, raw_log2FC, norm_log2FC, baseMean,
         raw_padj, norm_padj,
         eQTL_found, eQTL_count, mean_slope, negative_slope_count,
         eQTL_explains_expression)

write_csv(eqtl_results, "results/tables/eqtl_matches.csv")
cat("  Saved: results/tables/eqtl_matches.csv\n")

# Save final classification
write_csv(final_classification,
          "results/tables/final_gene_classification.csv")
cat("  Saved: results/tables/final_gene_classification.csv\n")

# =============================================================================
# STEP 7: Create summary visualization
# =============================================================================

cat("\nStep 7: Creating eQTL summary visualization...\n")

# Bar plot of eQTL findings
eqtl_summary_data <- eqtl_matches %>%
  count(eQTL_explains_expression) %>%
  mutate(
    percentage = 100 * n / sum(n),
    label = sprintf("%s\n(n=%d, %.0f%%)",
                    eQTL_explains_expression, n, percentage)
  )

p_eqtl <- ggplot(eqtl_summary_data,
                 aes(x = reorder(eQTL_explains_expression, -n),
                     y = n, fill = eQTL_explains_expression)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = n), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(
    values = c("eQTL supports low expression" = "#4DAF4A",
               "No eQTL found" = "#999999",
               "eQTL does not support low expression" = "#E41A1C")
  ) +
  labs(
    title = "eQTL Analysis of Chr21 DE Genes (<1.5 FC)",
    subtitle = paste0("Total genes analyzed: ", n_de_genes),
    x = "",
    y = "Number of Genes"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  )

ggsave("results/figures/eqtl_summary.pdf", p_eqtl,
       width = 8, height = 6)
cat("  Saved: results/figures/eqtl_summary.pdf\n")

# =============================================================================
# STEP 8: Final summary
# =============================================================================

cat("\n=== eQTL Analysis Summary ===\n")

cat(sprintf("\nChr21 genes analyzed: %d\n", nrow(chr21_categorized)))
cat(sprintf("DE genes with FC < 1.5: %d\n", n_de_genes))

if (n_de_genes > 0) {
  cat(sprintf("\neQTL findings:\n"))
  cat(sprintf("  Genes with eQTLs: %d (%.1f%%)\n",
              n_with_eqtl, 100 * n_with_eqtl / n_de_genes))
  cat(sprintf("  Low expression explained by eQTLs: %d\n", n_explained))
  cat(sprintf("  Low expression NOT explained: %d\n", n_unexplained))

  if (n_explained > 0) {
    cat("\n  Genes with eQTL support:\n")
    explained_genes <- eqtl_matches %>%
      filter(eQTL_explains_expression == "eQTL supports low expression") %>%
      select(Gene_name, eQTL_count, mean_slope)
    print(as.data.frame(explained_genes))
  }

  if (n_unexplained > 0) {
    cat("\n  Genes WITHOUT eQTL explanation:\n")
    unexplained_genes <- eqtl_matches %>%
      filter(eQTL_explains_expression != "eQTL supports low expression") %>%
      select(Gene_name, eQTL_found, eQTL_explains_expression)
    print(as.data.frame(unexplained_genes))
  }
}

# Final interpretation
cat("\n=== Final Interpretation ===\n")
total_genes <- nrow(chr21_categorized)
expected_dosage <- sum(chr21_categorized$category ==
                       "Expected Dosage (>=1.5 FC)")
pct_expected <- 100 * expected_dosage / total_genes

cat(sprintf("%.1f%% of chr21 genes (%d/%d) show expected 1.5-fold change.\n",
            pct_expected, expected_dosage, total_genes))

if (n_de_genes == 0) {
  cat("\nNo genes show significant deviation from expected dosage.\n")
  cat("This provides strong evidence AGAINST dosage compensation.\n")
} else if (n_explained >= n_de_genes * 0.5) {
  cat(sprintf("\n%d/%d genes with low expression are explained by eQTLs.\n",
              n_explained, n_de_genes))
  cat("Apparent 'dosage compensation' is due to genetic variation,\n")
  cat("NOT active compensation mechanisms.\n")
} else {
  cat(sprintf("\n%d genes show unexplained low expression.\n", n_unexplained))
  cat("These may represent true dosage compensation or\n")
  cat("could be explained by unknown eQTLs or other mechanisms.\n")
}

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/tables/eqtl_session_info.txt")

cat("\n=== Analysis Complete ===\n")
cat("All results saved to results/ directory.\n\n")
