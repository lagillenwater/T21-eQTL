# 02_categorize_genes.R
#
# Purpose: Categorize chromosome 21 genes based on fold change and expression
#          Reproduces the categorization logic from Hunter et al. Panel D
#
# Inputs:
#   - results/tables/deseq2_chr21_genes_both_analyses.csv
#   - data/processed/gene_annotations.csv
#   - data/processed/blacklisted_genes.csv (from process_blacklist.R)
#
# Outputs:
#   - results/tables/chr21_genes_categorized.csv
#   - results/tables/category_summary.csv
#
# Date: 2025-11-11

# Load required libraries
library(tidyverse)

# Set seed for reproducibility
set.seed(42)

cat("=== T21-eQTL Analysis: Chr21 Gene Categorization ===\n\n")

# =============================================================================
# STEP 1: Load results
# =============================================================================

cat("Step 1: Loading DESeq2 results...\n")

# Load the combined results file with both raw and normalized analyses
if (!file.exists("results/tables/deseq2_chr21_genes_both_analyses.csv")) {
  stop("Chr21 results not found. Run 01_deseq2_analysis.R first.")
}

chr21_results <- read_csv("results/tables/deseq2_chr21_genes_both_analyses.csv",
                          show_col_types = FALSE)

cat(sprintf("  Loaded %d chr21 genes (all types)\n", nrow(chr21_results)))

# Filter to protein-coding genes only (matching paper methodology)
# Paper reported 262 genes; they likely used protein-coding only
chr21_results <- chr21_results %>%
  filter(Gene_type == "protein_coding")

cat(sprintf("  Filtered to %d protein-coding genes\n", nrow(chr21_results)))

# =============================================================================
# STEP 2: Define categorization parameters
# =============================================================================

cat("\nStep 2: Defining categorization parameters...\n")

# Fold change thresholds
# Before ploidy normalization: chr21 genes expected at FC = 1.5 (log2FC = 0.585)
# After ploidy normalization: chr21 genes expected at FC = 1.0 (log2FC = 0)
FC_THRESHOLD <- log2(1.5)  # 0.585

# Expression threshold: use 2nd quintile (20th percentile) of baseMean
# Genes below this threshold have unreliable fold change estimates
# Paper: "set an expression level cutoff at the second quintile of the baseMean"
expr_threshold <- quantile(chr21_results$baseMean, 0.2, na.rm = TRUE)

# "No Reads" threshold: genes with essentially no expression
# These are genes where counts are too low to analyze
NO_READS_THRESHOLD <- 1  # baseMean < 1 means essentially no reads

# Significance threshold for ploidy-normalized analysis
PADJ_THRESHOLD <- 0.05

cat(sprintf("  FC threshold: log2FC >= %.3f (1.5-fold)\n", FC_THRESHOLD))
cat(sprintf("  No reads threshold: baseMean < %.1f\n", NO_READS_THRESHOLD))
cat(sprintf("  Low expression threshold: baseMean < %.2f (2nd quintile)\n",
            expr_threshold))
cat(sprintf("  Significance threshold: padj < %.2f\n", PADJ_THRESHOLD))

# =============================================================================
# STEP 3: Identify high genomic repeat genes
# =============================================================================

cat("\nStep 3: Identifying high genomic repeat genes...\n")

# Load ENCODE blacklist-derived gene list (from process_blacklist.R)
# These genes overlap with ENCODE blacklist regions (high signal/low mappability)
blacklist_file <- "data/processed/blacklisted_genes.csv"

if (file.exists(blacklist_file)) {
  blacklist_genes <- read_csv(blacklist_file, show_col_types = FALSE)
  blacklist_gene_names <- blacklist_genes$Gene_name
  cat(sprintf("  Loaded %d genes from ENCODE blacklist overlap analysis\n",
              length(blacklist_gene_names)))
} else {
  cat("  WARNING: Blacklist file not found. Run scripts/process_blacklist.R first.\n")
  cat("  Using fallback list of known repeat-rich genes.\n")
  blacklist_gene_names <- character(0)
}

# Additional known repeat-rich genes on chr21 (from literature/RepeatMasker)
# These include genes with paralogs or high repeat content
KNOWN_REPEAT_GENES <- c(
  # Ribosomal genes with paralogs
  "RPS6KB1", "RPS27", "RPS27L", "RPS27P",
  # Interferon receptor genes (some have paralogs)
  "IFNAR1", "IFNAR2",
  # Other known repeat-rich genes or genes with paralogs
  "TPTE", "BAGE", "DAB1"
)

# Combine ENCODE blacklist genes with known repeat-rich genes
HIGH_REPEAT_GENES <- unique(c(blacklist_gene_names, KNOWN_REPEAT_GENES))

chr21_results <- chr21_results %>%
  mutate(high_genomic_repeats = Gene_name %in% HIGH_REPEAT_GENES)

n_high_repeat <- sum(chr21_results$high_genomic_repeats, na.rm = TRUE)
cat(sprintf("  Total high repeat genes identified: %d\n", n_high_repeat))

if (n_high_repeat > 0) {
  cat("  Genes flagged:\n")
  high_repeat_list <- chr21_results %>%
    filter(high_genomic_repeats) %>%
    pull(Gene_name)
  cat("   ", paste(high_repeat_list, collapse = ", "), "\n")

  # Show source breakdown
  n_from_blacklist <- sum(high_repeat_list %in% blacklist_gene_names)
  n_from_known <- sum(high_repeat_list %in% KNOWN_REPEAT_GENES)
  cat(sprintf("  Source: %d from ENCODE blacklist, %d from known repeat genes\n",
              n_from_blacklist, n_from_known))
}

# =============================================================================
# STEP 4: Categorize genes (matching Panel D from Hunter et al.)
# =============================================================================

cat("\nStep 4: Categorizing chr21 genes (Panel D workflow)...\n")

#' Categorize a chr21 gene following Panel D workflow
#'
#' Panel D Flow (Hunter et al. 2023):
#' 1. First split: >= 1.5 FC vs < 1.5 FC (using raw analysis)
#' 2. For < 1.5 FC genes, categorize in order:
#'    a. No Reads - genes with essentially zero expression
#'    b. High Genomic Repeats - genes with multimapping issues
#'    c. Low Gene Expression - genes below 2nd quintile baseMean
#'    d. Not DE after ploidy normalization - genes where low FC is explained
#'       by proper ploidy accounting (brown branch in Panel D)
#'    e. DE after ploidy normalization - true candidates for dosage compensation
#'
#' @param raw_log2fc Raw log2 fold change (no ploidy normalization)
#' @param norm_padj Adjusted p-value from ploidy-normalized analysis
#' @param basemean Mean expression level
#' @param high_repeat Logical, whether gene has high genomic repeats
#' @param fc_thresh Fold change threshold (log2 scale, default log2(1.5))
#' @param no_reads_thresh Threshold for "no reads" category
#' @param expr_thresh Expression threshold for "low expression"
#' @param padj_thresh Significance threshold
#'
#' @return Category string
categorize_chr21_gene <- function(raw_log2fc, norm_padj, basemean, high_repeat,
                                   fc_thresh = log2(1.5),
                                   no_reads_thresh = 1,
                                   expr_thresh = 10,
                                   padj_thresh = 0.05) {

  # Handle NA values in key fields
  if (is.na(basemean)) {
    return("Insufficient Data")
  }

  # FIRST SPLIT: >= 1.5 FC vs < 1.5 FC (using raw/uncorrected analysis)
  # Genes at expected dosage show FC >= 1.5 (3 copies / 2 copies)
  if (!is.na(raw_log2fc) && raw_log2fc >= fc_thresh) {
    return("Expected Dosage (>=1.5 FC)")
  }

  # FOR GENES < 1.5 FC: Categorize by reason for low fold change


  # Category: No Reads
  # Genes with essentially no expression (cannot reliably measure FC)
  if (basemean < no_reads_thresh) {
    return("No Reads")
  }

  # Category: High Genomic Repeats
  # Genes with multimapping issues causing artifactually low counts
  if (high_repeat) {
    return("High Genomic Repeats")
  }

  # Category: Low Gene Expression
  # Genes below expression threshold have high variance in FC estimates
  if (basemean < expr_thresh) {
    return("Low Gene Expression")
  }

  # PLOIDY-NORMALIZED ANALYSIS: Is the gene DE after accounting for ploidy?
  # If padj >= 0.05: Not significantly different from expected (brown branch)
  # If padj < 0.05: Truly differentially expressed (candidates for compensation)

  if (is.na(norm_padj) || norm_padj >= padj_thresh) {
    # Not significantly DE after ploidy normalization
    # The low FC is explained by proper ploidy accounting
    return("Expression/Ploidy Not Differentially Expressed")
  } else {
    # Significantly DE even after ploidy normalization
    # These are candidates for dosage compensation or eQTL effects
    return("Expression/Ploidy Differentially Expressed")
  }
}

# Apply categorization
chr21_categorized <- chr21_results %>%
  mutate(
    category = mapply(
      categorize_chr21_gene,
      raw_log2fc = raw_log2FC,
      norm_padj = norm_padj,
      basemean = baseMean,
      high_repeat = high_genomic_repeats,
      fc_thresh = FC_THRESHOLD,
      no_reads_thresh = NO_READS_THRESHOLD,
      expr_thresh = expr_threshold,
      padj_thresh = PADJ_THRESHOLD
    ),
    # Calculate actual fold changes for reference
    raw_fold_change = 2^raw_log2FC,
    norm_fold_change = 2^norm_log2FC
  )

# =============================================================================
# STEP 5: Generate category summary
# =============================================================================

cat("\nStep 5: Generating category summary...\n")

category_summary <- chr21_categorized %>%
  count(category) %>%
  arrange(desc(n)) %>%
  mutate(percentage = 100 * n / sum(n))

cat("\n  Category breakdown:\n")
print(as.data.frame(category_summary))

# Compare with paper results (for reference)
cat("\n  Expected from Hunter et al. (2023) Panel D:\n")
cat("    Expected Dosage (>=1.5 FC): 101 genes\n")
cat("    <1.5 FC: 151 genes, broken down as:\n")
cat("      - No Reads: 10 genes\n")
cat("      - High Genomic Repeats: 9 genes\n")
cat("      - Low Gene Expression: 58 genes\n")
cat("      - Expression/Ploidy Not DE: 79 genes (brown branch)\n")
cat("      - Expression/Ploidy DE: 5 genes (candidates for compensation)\n")
cat("    Total: 262 genes (after filtering)\n")

# =============================================================================
# STEP 6: Detailed analysis of DE genes after ploidy normalization
# =============================================================================

cat("\nStep 6: Analyzing Expression/Ploidy Differentially Expressed genes...\n")
cat("  (These are genes with <1.5 FC that remain DE after ploidy normalization)\n")

de_genes_low_fc <- chr21_categorized %>%
  filter(category == "Expression/Ploidy Differentially Expressed") %>%
  arrange(norm_padj)

n_de_low_fc <- nrow(de_genes_low_fc)
cat(sprintf("  Found %d genes\n", n_de_low_fc))

if (n_de_low_fc > 0) {
  cat("\n  These genes are candidates for dosage compensation or eQTL effects:\n")
  de_summary <- de_genes_low_fc %>%
    select(Gene_name, raw_log2FC, raw_fold_change, norm_log2FC,
           norm_fold_change, baseMean, raw_padj, norm_padj) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  print(as.data.frame(de_summary))

  cat("\n  These genes will be cross-referenced with eQTL databases\n")
  cat("  to determine if expression level is supported by canonical eQTLs.\n")
  cat("  See 05_eqtl_analysis.R\n")
} else {
  cat("\n  No genes remain differentially expressed after ploidy normalization\n")
  cat("  This suggests genes with <1.5 FC are explained by:\n")
  cat("    - Technical artifacts (no reads, repeats, low expression)\n")
  cat("    - Proper accounting for ploidy differences\n")
}

# =============================================================================
# STEP 7: Save results
# =============================================================================

cat("\nStep 7: Saving categorized results...\n")

# Save full categorization
chr21_categorized %>%
  select(EnsemblID, Gene_name, Chr, Gene_type,
         baseMean, raw_log2FC, raw_fold_change, raw_pvalue, raw_padj,
         norm_log2FC, norm_fold_change, norm_pvalue, norm_padj,
         high_genomic_repeats, category) %>%
  arrange(category, raw_padj) %>%
  write_csv("results/tables/chr21_genes_categorized.csv")

cat("  Saved: results/tables/chr21_genes_categorized.csv\n")

# Save category summary
category_summary %>%
  write_csv("results/tables/category_summary.csv")

cat("  Saved: results/tables/category_summary.csv\n")

# Save list of DE genes for eQTL analysis
if (n_de_low_fc > 0) {
  de_genes_low_fc %>%
    select(EnsemblID, Gene_name, raw_log2FC, norm_log2FC, baseMean,
           raw_padj, norm_padj) %>%
    write_csv("results/tables/chr21_de_genes_for_eqtl.csv")
  cat("  Saved: results/tables/chr21_de_genes_for_eqtl.csv\n")
}

# =============================================================================
# STEP 8: Create category flow diagram data (Panel D structure)
# =============================================================================

cat("\nStep 8: Preparing data for alluvial plot (Panel D)...\n")

# Create flow data matching Panel D structure:
# Level 1: Chr21 Genes (total)
# Level 2: >=1.5 FC vs <1.5 FC
# Level 3 (for <1.5 FC): No Reads, High Genomic Repeats, Low Expression,
#                        Ploidy Not DE, Ploidy DE
# Level 4 (for Ploidy DE): eQTL supports vs doesn't support (added in eQTL analysis)

flow_data <- chr21_categorized %>%
  mutate(
    # First split: FC threshold (using raw log2FC)
    fc_group = ifelse(raw_log2FC >= log2(1.5),
                      ">=1.5 FC", "<1.5 FC"),
    # Second split: detailed category for <1.5 FC genes
    detailed_category = case_when(
      fc_group == ">=1.5 FC" ~ "Expected Dosage",
      category == "No Reads" ~ "No Reads",
      category == "High Genomic Repeats" ~ "High Genomic Repeats",
      category == "Low Gene Expression" ~ "Low Gene Expression",
      category == "Expression/Ploidy Not Differentially Expressed" ~
        "Expression/Ploidy Not DE",
      category == "Expression/Ploidy Differentially Expressed" ~
        "Expression/Ploidy DE",
      TRUE ~ "Other"
    )
  ) %>%
  count(fc_group, detailed_category)

cat("\n  Flow diagram structure (Panel D):\n")
print(as.data.frame(flow_data))

flow_data %>%
  write_csv("results/tables/alluvial_flow_data.csv")
cat("  Saved: results/tables/alluvial_flow_data.csv\n")

# =============================================================================
# STEP 9: Summary
# =============================================================================

cat("\n=== Categorization Summary ===\n")
cat(sprintf("Total chr21 genes analyzed: %d\n", nrow(chr21_categorized)))
cat(sprintf("\nCategory counts:\n"))

for (i in 1:nrow(category_summary)) {
  cat(sprintf("  %s: %d (%.1f%%)\n",
              category_summary$category[i],
              category_summary$n[i],
              category_summary$percentage[i]))
}

cat(sprintf("\nGenes requiring eQTL analysis: %d\n", n_de_low_fc))

# Interpretation
cat("\n=== Interpretation ===\n")

n_expected <- category_summary %>%
  filter(category == "Expected Dosage (>=1.5 FC)") %>%
  pull(n)

n_ploidy_de <- category_summary %>%
  filter(category == "Expression/Ploidy Differentially Expressed") %>%
  pull(n)

if (length(n_expected) == 0) n_expected <- 0
if (length(n_ploidy_de) == 0) n_ploidy_de <- 0

cat(sprintf("Genes at expected dosage (>=1.5 FC): %d\n", n_expected))
cat(sprintf("Genes DE after ploidy normalization: %d\n", n_ploidy_de))

if (n_ploidy_de <= 10) {
  cat("\nFew genes remain DE after proper ploidy accounting.\n")
  cat("This is consistent with NO widespread dosage compensation.\n")
  cat("Genes with <1.5 FC are explained by:\n")
  cat("  - Technical artifacts (no reads, repeats, low expression)\n")
  cat("  - Statistical noise resolved by ploidy normalization\n")
  if (n_ploidy_de > 0) {
    cat(sprintf("  - %d genes require eQTL analysis to explain\n", n_ploidy_de))
  }
} else {
  cat("\nMultiple genes remain DE after ploidy normalization.\n")
  cat("eQTL analysis needed to determine if these represent:\n")
  cat("  - True dosage compensation\n")
  cat("  - Allelic variation (eQTL effects)\n")
}

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/tables/categorization_session_info.txt")
cat("\nSaved session info\n")

cat("\n=== Categorization Complete ===\n")
cat("Next steps:\n")
cat("  - Run 03_volcano_plot.R for visualization\n")
cat("  - Run 04_alluvial_plot.R to create Panel D\n")
cat("  - Run 05_eqtl_analysis.R for eQTL cross-reference\n\n")
