# 07_sankeymatic_export.R
#
# Purpose: Generate SankeyMATIC input file for chr21 gene classification diagram
#          Creates flow data for visualization at https://sankeymatic.com/build/
#          Includes both < 1.5 FC analysis (potential compensation) and
#          >= 1.5 FC analysis (expected dosage genes significant after ploidy)
#
# Inputs:
#   - results/tables/final_gene_classification_all.csv (with eQTL data for all genes)
#   - results/tables/expected_dosage_eqtl_matches.csv (eQTL analysis for >= 1.5 FC genes)
#   - data/processed/blacklisted_genes.csv
#
# Outputs:
#   - results/sankeymatic_input.txt
#   - results/sankeymatic_input_extended.txt (includes >= 1.5 FC breakdown)
#
# Author: Claude Code
# Date: 2025-12-17

library(tidyverse)

set.seed(42)

cat("=== T21-eQTL Analysis: SankeyMATIC Export ===\n\n")

# =============================================================================
# STEP 1: Load data
# =============================================================================

cat("Step 1: Loading data...\n")

if (!file.exists("results/tables/final_gene_classification_all.csv")) {
  stop("Gene classification file not found. Run 05_eqtl_analysis.R first.")
}

chr21_full <- read_csv("results/tables/final_gene_classification_all.csv",
                       show_col_types = FALSE)

cat(sprintf("  Loaded %d chr21 genes (protein_coding + lncRNA)\n", nrow(chr21_full)))

# Load expected dosage eQTL results if available
expected_dosage_file <- "results/tables/expected_dosage_eqtl_matches.csv"
if (file.exists(expected_dosage_file)) {
  expected_dosage_eqtl <- read_csv(expected_dosage_file, show_col_types = FALSE)
  cat(sprintf("  Loaded %d expected dosage genes with eQTL analysis\n",
              nrow(expected_dosage_eqtl)))
  has_expected_dosage_analysis <- TRUE
} else {
  cat("  Expected dosage eQTL file not found - run 08_expected_dosage_eqtl.R\n")
  has_expected_dosage_analysis <- FALSE
}

# =============================================================================
# STEP 2: Load high repeat gene list
# =============================================================================

cat("\nStep 2: Loading high repeat gene annotations...\n")

blacklist_file <- "data/processed/blacklisted_genes.csv"
if (file.exists(blacklist_file)) {
  blacklist_genes <- read_csv(blacklist_file, show_col_types = FALSE)
  blacklist_gene_names <- blacklist_genes$Gene_name
  cat(sprintf("  Loaded %d blacklisted genes\n", length(blacklist_gene_names)))
} else {
  blacklist_gene_names <- character(0)
  cat("  No blacklist file found\n")
}

# Known repeat-rich genes from literature
KNOWN_REPEAT_GENES <- c(
  "RPS6KB1", "RPS27", "RPS27L", "RPS27P",
  "IFNAR1", "IFNAR2",
  "TPTE", "BAGE", "DAB1"
)

HIGH_REPEAT_GENES <- unique(c(blacklist_gene_names, KNOWN_REPEAT_GENES))
cat(sprintf("  Total high repeat genes: %d\n", length(HIGH_REPEAT_GENES)))

# Add high repeat flag
chr21_full <- chr21_full %>%
  mutate(high_genomic_repeats = Gene_name %in% HIGH_REPEAT_GENES)

# =============================================================================
# STEP 3: Calculate thresholds and categorize genes
# =============================================================================

cat("\nStep 3: Categorizing genes...\n")

# Low expression threshold: 2nd quintile of baseMean
basemean_threshold <- quantile(chr21_full$baseMean, 0.2)
cat(sprintf("  Low expression threshold (20th percentile): %.2f\n", basemean_threshold))

# Split by raw FC (before ploidy correction)
# Expected dosage: raw FC >= 1.5
n_above_15 <- sum(chr21_full$raw_log2FC >= log2(1.5))
n_below_15 <- sum(chr21_full$raw_log2FC < log2(1.5))

cat(sprintf("  Genes >= 1.5 FC (Expected Dosage): %d\n", n_above_15))
cat(sprintf("  Genes < 1.5 FC: %d\n", n_below_15))

# =============================================================================
# STEP 4: Categorize < 1.5 FC genes
# =============================================================================

cat("\nStep 4: Categorizing < 1.5 FC genes...\n")

# Priority: High Repeats > Low Expression > Not DE > DE
below_15 <- chr21_full %>%
  filter(raw_log2FC < log2(1.5)) %>%
  mutate(
    category = case_when(
      high_genomic_repeats ~ "High Repeats",
      baseMean < basemean_threshold ~ "Low Expression",
      norm_padj >= 0.05 | is.na(norm_padj) ~ "Not DE",
      TRUE ~ "DE"
    )
  )

n_high_rep <- sum(below_15$category == "High Repeats")
n_low_expr <- sum(below_15$category == "Low Expression")
n_not_de <- sum(below_15$category == "Not DE")
n_de <- sum(below_15$category == "DE")

cat(sprintf("  < 1.5 FC breakdown:\n"))
cat(sprintf("    High Repeats: %d\n", n_high_rep))
cat(sprintf("    Low Expression: %d\n", n_low_expr))
cat(sprintf("    Not DE (padj >= 0.05): %d\n", n_not_de))
cat(sprintf("    DE (padj < 0.05): %d\n", n_de))

# eQTL classification for < 1.5 FC DE genes
de_genes <- below_15 %>% filter(category == "DE")
n_de_has_eqtl <- sum(de_genes$eQTL_found == TRUE, na.rm = TRUE)
n_de_no_eqtl <- sum(de_genes$eQTL_found == FALSE | is.na(de_genes$eQTL_found))

cat(sprintf("\n  < 1.5 FC DE genes eQTL status:\n"))
cat(sprintf("    Has eQTL: %d\n", n_de_has_eqtl))
cat(sprintf("    No eQTL: %d\n", n_de_no_eqtl))

# =============================================================================
# STEP 5: Categorize >= 1.5 FC genes (Extended Analysis)
# =============================================================================

cat("\nStep 5: Categorizing >= 1.5 FC genes (Extended Analysis)...\n")

above_15 <- chr21_full %>%
  filter(raw_log2FC >= log2(1.5))

# Split by significance after ploidy correction
n_above_15_not_sig <- sum(above_15$norm_padj >= 0.05 | is.na(above_15$norm_padj))
n_above_15_sig <- sum(above_15$norm_padj < 0.05, na.rm = TRUE)

cat(sprintf("  >= 1.5 FC breakdown:\n"))
cat(sprintf("    Not significant after ploidy (as expected): %d\n", n_above_15_not_sig))
cat(sprintf("    Still significant after ploidy (deviation): %d\n", n_above_15_sig))

# eQTL classification for significant >= 1.5 FC genes (simple: has eQTL vs no eQTL)
if (has_expected_dosage_analysis && n_above_15_sig > 0) {
  # Simple classification: has eQTL or not
  n_exp_has_eqtl <- sum(expected_dosage_eqtl$eQTL_found == TRUE, na.rm = TRUE)
  n_exp_no_eqtl <- sum(expected_dosage_eqtl$eQTL_found == FALSE |
                         is.na(expected_dosage_eqtl$eQTL_found))

  # Calculate genes not in eQTL analysis (e.g., lncRNAs not analyzed)
  n_exp_total_analyzed <- n_exp_has_eqtl + n_exp_no_eqtl
  n_exp_not_analyzed <- n_above_15_sig - n_exp_total_analyzed

  cat(sprintf("\n  >= 1.5 FC significant genes eQTL status:\n"))
  cat(sprintf("    Has eQTL: %d\n", n_exp_has_eqtl))
  cat(sprintf("    No eQTL: %d\n", n_exp_no_eqtl))
  if (n_exp_not_analyzed > 0) {
    cat(sprintf("    Not in eQTL analysis (lncRNAs, etc.): %d\n", n_exp_not_analyzed))
  }
} else {
  n_exp_has_eqtl <- 0
  n_exp_no_eqtl <- n_above_15_sig
  n_exp_not_analyzed <- 0
}

# =============================================================================
# STEP 6: Generate SankeyMATIC input (Original format)
# =============================================================================

cat("\nStep 6: Generating SankeyMATIC input (original format)...\n")

# Build the original Sankey flow text (for backwards compatibility)
sankey_lines_original <- c(
  sprintf("Chr21 Genes [%d] >=1.5 FC", n_above_15),
  sprintf("Chr21 Genes [%d] <1.5 FC", n_below_15),
  sprintf(">=1.5 FC [%d] Expected Dosage", n_above_15),
  sprintf("<1.5 FC [%d] High Repeats", n_high_rep),
  sprintf("<1.5 FC [%d] Low Expression", n_low_expr),
  sprintf("<1.5 FC [%d] Not DE", n_not_de),
  sprintf("<1.5 FC [%d] DE", n_de),
  sprintf("DE [%d] Has eQTL", n_de_has_eqtl),
  sprintf("DE [%d] No eQTL", n_de_no_eqtl)
)

# Write original format
output_file_original <- "results/sankeymatic_input.txt"
writeLines(sankey_lines_original, output_file_original)
cat(sprintf("  Saved: %s\n", output_file_original))

# =============================================================================
# STEP 7: Generate SankeyMATIC input (Extended format with >= 1.5 FC analysis)
# =============================================================================

cat("\nStep 7: Generating SankeyMATIC input (extended format)...\n")

# Build extended Sankey flow with both branches fully analyzed
sankey_lines_extended <- c(
  "// Chr21 Gene Expression Analysis - Extended",
  "// Generated by 07_sankeymatic_export.R",
  "",
  "// Main split by fold change",
  sprintf("Chr21 Genes [%d] >=1.5 FC", n_above_15),
  sprintf("Chr21 Genes [%d] <1.5 FC", n_below_15),
  "",
  "// >= 1.5 FC branch (Expected Dosage Analysis)",
  sprintf(">=1.5 FC [%d] Not Sig After Ploidy", n_above_15_not_sig),
  sprintf(">=1.5 FC [%d] Sig After Ploidy", n_above_15_sig),
  sprintf("Not Sig After Ploidy [%d] Expected Dosage", n_above_15_not_sig)
)

# Add eQTL breakdown for significant >= 1.5 FC genes
if (n_above_15_sig > 0) {
  sankey_lines_extended <- c(
    sankey_lines_extended,
    sprintf("Sig After Ploidy [%d] Has eQTL", n_exp_has_eqtl),
    sprintf("Sig After Ploidy [%d] No eQTL", n_exp_no_eqtl)
  )
  if (n_exp_not_analyzed > 0) {
    sankey_lines_extended <- c(
      sankey_lines_extended,
      sprintf("Sig After Ploidy [%d] Not Analyzed", n_exp_not_analyzed)
    )
  }
}

# Add < 1.5 FC branch
sankey_lines_extended <- c(
  sankey_lines_extended,
  "",
  "// < 1.5 FC branch (Potential Compensation Analysis)",
  sprintf("<1.5 FC [%d] High Repeats", n_high_rep),
  sprintf("<1.5 FC [%d] Low Expression", n_low_expr),
  sprintf("<1.5 FC [%d] Not DE", n_not_de),
  sprintf("<1.5 FC [%d] DE Low FC", n_de),
  sprintf("DE Low FC [%d] Has eQTL", n_de_has_eqtl),
  sprintf("DE Low FC [%d] No eQTL", n_de_no_eqtl)
)

# Add color suggestions
sankey_lines_extended <- c(
  sankey_lines_extended,
  "",
  "// Color suggestions:",
  "// :Chr21 Genes #808080",
  "// :>=1.5 FC #4393c3",
  "// :<1.5 FC #d6604d",
  "// :Not Sig After Ploidy #92c5de",
  "// :Sig After Ploidy #2166ac",
  "// :Expected Dosage #67a9cf",
  "// :Has eQTL #4DAF4A",
  "// :No eQTL #999999",
  "// :Not Analyzed #CCCCCC",
  "// :High Repeats #f4a582",
  "// :Low Expression #fddbc7",
  "// :Not DE #d1e5f0",
  "// :DE Low FC #b2182b"
)

# Write extended format
output_file_extended <- "results/sankeymatic_input_extended.txt"
writeLines(sankey_lines_extended, output_file_extended)
cat(sprintf("  Saved: %s\n", output_file_extended))

# =============================================================================
# STEP 8: Generate clean version for direct paste
# =============================================================================

cat("\nStep 8: Generating clean paste-ready version...\n")

# Clean version without comments for direct paste
sankey_clean <- c(
  sprintf("Chr21 Genes [%d] >=1.5 FC", n_above_15),
  sprintf("Chr21 Genes [%d] <1.5 FC", n_below_15),
  sprintf(">=1.5 FC [%d] Not Sig After Ploidy", n_above_15_not_sig),
  sprintf(">=1.5 FC [%d] Sig After Ploidy", n_above_15_sig),
  sprintf("Not Sig After Ploidy [%d] Expected Dosage", n_above_15_not_sig),
  sprintf("Sig After Ploidy [%d] Has eQTL", n_exp_has_eqtl),
  sprintf("Sig After Ploidy [%d] No eQTL", n_exp_no_eqtl)
)

# Add "Not Analyzed" if there are genes missing from eQTL analysis
if (n_exp_not_analyzed > 0) {
  sankey_clean <- c(
    sankey_clean,
    sprintf("Sig After Ploidy [%d] Not Analyzed", n_exp_not_analyzed)
  )
}

# Add < 1.5 FC branch
sankey_clean <- c(
  sankey_clean,
  sprintf("<1.5 FC [%d] High Repeats", n_high_rep),
  sprintf("<1.5 FC [%d] Low Expression", n_low_expr),
  sprintf("<1.5 FC [%d] Not DE", n_not_de),
  sprintf("<1.5 FC [%d] DE Low FC", n_de),
  sprintf("DE Low FC [%d] Has eQTL", n_de_has_eqtl),
  sprintf("DE Low FC [%d] No eQTL", n_de_no_eqtl)
)

output_file_clean <- "results/sankeymatic_input_clean.txt"
writeLines(sankey_clean, output_file_clean)
cat(sprintf("  Saved: %s\n", output_file_clean))

# =============================================================================
# STEP 9: Summary statistics
# =============================================================================

cat("\n=== SUMMARY ===\n")
cat(sprintf("Total chr21 genes (protein_coding + lncRNA): %d\n", nrow(chr21_full)))

cat("\n--- >= 1.5 FC Branch (Expected Dosage) ---\n")
cat(sprintf("  Total: %d (%.1f%%)\n", n_above_15, 100 * n_above_15 / nrow(chr21_full)))
cat(sprintf("  Not significant after ploidy (truly expected): %d (%.1f%%)\n",
            n_above_15_not_sig, 100 * n_above_15_not_sig / n_above_15))
cat(sprintf("  Significant after ploidy (deviates from expected): %d (%.1f%%)\n",
            n_above_15_sig, 100 * n_above_15_sig / n_above_15))

if (n_above_15_sig > 0) {
  cat(sprintf("    - Has eQTL: %d (%.1f%%)\n",
              n_exp_has_eqtl, 100 * n_exp_has_eqtl / n_above_15_sig))
  cat(sprintf("    - No eQTL: %d (%.1f%%)\n",
              n_exp_no_eqtl, 100 * n_exp_no_eqtl / n_above_15_sig))
  if (n_exp_not_analyzed > 0) {
    cat(sprintf("    - Not analyzed (lncRNAs, etc.): %d (%.1f%%)\n",
                n_exp_not_analyzed, 100 * n_exp_not_analyzed / n_above_15_sig))
  }
}

cat("\n--- < 1.5 FC Branch (Potential Compensation) ---\n")
cat(sprintf("  Total: %d (%.1f%%)\n", n_below_15, 100 * n_below_15 / nrow(chr21_full)))
cat(sprintf("  Technical/filtering:\n"))
cat(sprintf("    - High Repeats: %d\n", n_high_rep))
cat(sprintf("    - Low Expression: %d\n", n_low_expr))
cat(sprintf("  Not significant after ploidy: %d\n", n_not_de))
cat(sprintf("  Significant after ploidy (potential compensation): %d\n", n_de))

if (n_de > 0) {
  cat(sprintf("    - Has eQTL: %d (%.1f%%)\n", n_de_has_eqtl, 100 * n_de_has_eqtl / n_de))
  cat(sprintf("    - No eQTL: %d (%.1f%%)\n", n_de_no_eqtl, 100 * n_de_no_eqtl / n_de))
}

cat("\n=== Instructions ===\n")
cat("1. Go to https://sankeymatic.com/build/\n")
cat("2. Choose which file to use:\n")
cat("   - sankeymatic_input.txt: Original format (< 1.5 FC focus only)\n")
cat("   - sankeymatic_input_clean.txt: Extended format (both branches)\n")
cat("   - sankeymatic_input_extended.txt: With comments and color suggestions\n")
cat("3. Paste the contents and customize colors\n")
cat("4. Export as PNG/SVG\n")

cat("\n=== SankeyMATIC Input (Clean Extended) ===\n")
cat(paste(sankey_clean, collapse = "\n"), "\n")

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/sankeymatic_session_info.txt")

cat("\n=== SankeyMATIC Export Complete ===\n")
