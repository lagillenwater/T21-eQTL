# 12_chr21_lane_assignment.R
#
# Purpose: Assign every chr21 gene to a (significance lane, eQTL lane) pair
#          for the comprehensive Sankey/alluvial figure. Mirrors the paper's
#          framing (Hunter et al. 2023, BMC Biology 21:228): no raw fold-
#          change buffer, padj < 0.01 after ploidy normalization, and a
#          locus-level eQTL scan asking whether ANY cis variant in GTEx
#          whole blood matches the observed deviation direction. We also
#          require within-T21 reproducibility for a variant to "count" as
#          supporting - leveraging the n = 300+ cohort that the paper's
#          family-of-4 study could not.
#
#          Sig lane: DE_low (norm_padj < 0.01 AND norm_log2FC < 0)
#                    DE_high (norm_padj < 0.01 AND norm_log2FC > 0)
#                    Not_DE_or_NA (otherwise)
#          eQTL lane: explained        (>= 1 cis variant with both:
#                                       sign(gtex_slope) == sign(norm_log2FC)
#                                       AND sign(t21_slope) == sign(gtex_slope))
#                     unexplained      (cis variants exist but none satisfy)
#                     no_GTEx_data     (gene not in GTEx whole-blood
#                                       signif_pairs at all)
#
# Inputs:
#   - results/tables/deseq2_chr21_combined.csv
#   - results/tables/deseq2_all_genes_ploidy_normalized.csv  (cohort-noise SD)
#   - results/tables/t21_dosage_per_variant.csv
#   - data/processed/eqtl_target_variants.csv
#   - data/processed/blacklisted_genes.csv  (high-repeat flag)
#
# Outputs:
#   - results/tables/chr21_lane_assignments.csv
#   - results/tables/chr21_lane_summary.csv
#   - results/tables/chr21_lane_assignment_session_info.txt
#
# Author: Claude Code
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

cat("=== T21-eQTL: Per-Gene Lane Assignment (DE x eQTL) ===\n\n")

# =============================================================================
# Constants - matched to paper conventions
# =============================================================================

ALPHA               <- 0.01    # paper: padj < .01 (Fig. 2B, 3B)
ALPHA_REPRO         <- 0.05    # within-T21 nominal p for a cis variant to
                               # count as a reproducible supporter
MAGNITUDE_THRESHOLD <- 1.0     # |norm_log2FC| / cohort_sd >= this to be
                               # eligible for eQTL-based explanation. Genes
                               # below this threshold are flagged
                               # "below_cohort_noise" - their statistically
                               # significant deviations are within typical
                               # cohort variation, so no eQTL explanation
                               # is sought. This filter is also applied
                               # upstream in script 09 (gene selection)
                               # so failing genes have no genotype data.
LOW_EXPR_QUANT      <- 0.20    # paper: "second quintile of baseMean"
RESTRICT_TO_PROTEIN_CODING <- TRUE   # restrict chr21 set + cohort-noise
                                     # reference to protein-coding genes

KNOWN_REPEAT_GENES <- c("RPS6KB1", "RPS27", "RPS27L", "RPS27P",
                        "IFNAR1", "IFNAR2", "TPTE", "BAGE", "DAB1")

# =============================================================================
# STEP 1: Load inputs
# =============================================================================

cat("Step 1: Loading inputs...\n")

deseq <- fread("results/tables/deseq2_chr21_combined.csv")
deseq[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
if (RESTRICT_TO_PROTEIN_CODING) {
  n_before <- nrow(deseq)
  deseq <- deseq[Gene_type == "protein_coding"]
  cat(sprintf("  chr21 DESeq2 rows (protein_coding only): %d (was %d)\n",
              nrow(deseq), n_before))
} else {
  cat(sprintf("  chr21 DESeq2 rows: %d\n", nrow(deseq)))
}

per_var <- fread("results/tables/t21_dosage_per_variant.csv")
cat(sprintf("  Per-variant within-T21 fits: %d\n", nrow(per_var)))

targets <- fread("data/processed/eqtl_target_variants.csv")
cat(sprintf("  GTEx target (variant, gene) pairs: %d (%d unique genes)\n",
            nrow(targets), uniqueN(targets$ensembl_stable)))

blacklist_path <- "data/processed/blacklisted_genes.csv"
if (file.exists(blacklist_path)) {
  blacklist_genes <- fread(blacklist_path)$Gene_name
} else {
  blacklist_genes <- character(0)
  cat("  (no blacklisted_genes.csv found - using built-in repeat list only)\n")
}
high_repeat_genes <- unique(c(blacklist_genes, KNOWN_REPEAT_GENES))
cat(sprintf("  High-repeat genes flagged: %d\n", length(high_repeat_genes)))

# =============================================================================
# STEP 2: Per-gene locus-level aggregation of cis variants
# =============================================================================

cat("\nStep 2: Aggregating cis variants per gene...\n")

# Bring observed deviation sign onto each per-variant row
per_var <- merge(per_var,
                 deseq[, .(ensembl_stable, norm_log2FC_obs = norm_log2FC)],
                 by = "ensembl_stable", all.x = TRUE)

# Two flags per variant:
#   dir_match             - GTEx slope sign == observed deviation sign
#   supportive_with_repro - dir_match AND within-T21 slope reproduces GTEx
#                           in BOTH sign AND nominal significance
#                           (t21_p < ALPHA_REPRO). Sign-only matching
#                           lets noise through and produces unconvincing
#                           boxplot panels (e.g., APP at t21_p ~ 0.07).
per_var[, dir_match := !is.na(gtex_slope) & !is.na(norm_log2FC_obs) &
                       sign(gtex_slope) == sign(norm_log2FC_obs)]
per_var[, supportive_with_repro :=
          dir_match &
          !is.na(t21_slope) & !is.na(t21_p) &
          sign(t21_slope) == sign(gtex_slope) &
          t21_p < ALPHA_REPRO]

# Per-gene aggregation
locus <- per_var[, .(
  n_cis_total       = uniqueN(variant_id),
  n_dir_match       = uniqueN(variant_id[dir_match == TRUE]),
  n_supp_with_repro = uniqueN(variant_id[supportive_with_repro == TRUE])
), by = ensembl_stable]

# Strongest supportive-with-reproducibility variant per gene (smallest
# WITHIN-T21 p, since this variant is what the boxplot displays - we want
# the most visually convincing within-T21 example, not the most
# established GTEx eQTL). For unexplained genes, fall back to the
# strongest dir_match variant ranked by GTEx p (tells the boxplot script
# which variant to display as the "best try" that failed within T21).
pick_strongest <- function(dt, sort_col, pval_col, variant_col,
                           pval_out, variant_out) {
  setorderv(dt, sort_col)
  out <- dt[, list(.SD[[pval_col]][1L], .SD[[variant_col]][1L]),
            by = ensembl_stable,
            .SDcols = c(pval_col, variant_col)]
  setnames(out, c("V1", "V2"), c(pval_out, variant_out))
  out
}

strongest_repro <- pick_strongest(
  per_var[supportive_with_repro == TRUE], "t21_p", "t21_p", "variant_id",
  "strongest_supp_pval", "strongest_supp_variant")

strongest_dir <- pick_strongest(
  per_var[dir_match == TRUE], "gtex_pval", "gtex_pval", "variant_id",
  "strongest_dir_pval", "strongest_dir_variant")

# Fallback for unexplained genes with NO direction-matching variant: the
# strongest cis variant overall (smallest GTEx pval, regardless of direction).
# Gives the boxplot script something to display for every gene that has any
# cis variant tested.
strongest_overall <- pick_strongest(
  per_var, "gtex_pval", "gtex_pval", "variant_id",
  "strongest_overall_pval", "strongest_overall_variant")

locus <- merge(locus, strongest_repro,   by = "ensembl_stable", all.x = TRUE)
locus <- merge(locus, strongest_dir,     by = "ensembl_stable", all.x = TRUE)
locus <- merge(locus, strongest_overall, by = "ensembl_stable", all.x = TRUE)

cat(sprintf("  Genes with any cis variant tested in T21: %d\n", nrow(locus)))

# =============================================================================
# STEP 3: Build the per-gene table and assign lanes
# =============================================================================

cat("\nStep 3: Assigning lanes...\n")

m <- merge(deseq[, .(ensembl_stable, EnsemblID, Gene_name, Gene_type, Chr,
                     baseMean, raw_log2FC, raw_padj,
                     norm_log2FC, norm_padj)],
           locus,
           by = "ensembl_stable", all.x = TRUE)

# Genes with no GTEx whole-blood signif eQTL at all -> 0 cis variants tested
m[is.na(n_cis_total),       n_cis_total       := 0L]
m[is.na(n_dir_match),       n_dir_match       := 0L]
m[is.na(n_supp_with_repro), n_supp_with_repro := 0L]
m[, raw_FC := 2^raw_log2FC]

# Paper's filters as flags (categorize, do not remove)
basemean_threshold <- quantile(m$baseMean, LOW_EXPR_QUANT, na.rm = TRUE)
cat(sprintf("  baseMean cutoff at %.0fth percentile: %.2f\n",
            100 * LOW_EXPR_QUANT, basemean_threshold))
m[, low_expr   := !is.na(baseMean) & baseMean < basemean_threshold]
m[, high_repeat := Gene_name %in% high_repeat_genes]

# Deviation magnitude vs cohort noise (computed BEFORE lane assignment so
# the magnitude filter can gate the eQTL lane downstream).
all_genes_lfc <- fread(
  "results/tables/deseq2_all_genes_ploidy_normalized.csv")
if (RESTRICT_TO_PROTEIN_CODING) {
  all_genes_lfc <- all_genes_lfc[Gene_type == "protein_coding"]
}
non_chr21_sd <- sd(
  all_genes_lfc[Chr != "chr21" & !is.na(log2FoldChange), log2FoldChange])
cat(sprintf("  cohort-noise SD (non-chr21%s): %.3f\n",
            if (RESTRICT_TO_PROTEIN_CODING) ", protein-coding" else "",
            non_chr21_sd))
m[, deviation_magnitude := abs(norm_log2FC)]
m[, deviation_vs_cohort_sd := deviation_magnitude / non_chr21_sd]
m[, passes_magnitude_filter := !is.na(deviation_vs_cohort_sd) &
                               deviation_vs_cohort_sd >= MAGNITUDE_THRESHOLD]

# Significance lane - cohort-noise filter is the FIRST split. Genes whose
# ploidy-corrected deviation is within typical cohort variation are
# categorized as Expected_dosage regardless of raw FC sign or padj
# (statistical significance at large n is not the same as biological
# departure from ploidy expectation). Surviving genes get the standard
# categorization.
m[, sig_lane := fcase(
  passes_magnitude_filter == FALSE,                       "Expected_dosage",
  high_repeat == TRUE,                                    "High_repeats",
  low_expr == TRUE,                                       "Low_expression",
  !is.na(norm_padj) & norm_padj < ALPHA & norm_log2FC < 0, "DE_low",
  !is.na(norm_padj) & norm_padj < ALPHA & norm_log2FC > 0, "DE_high",
  default                                               = "Not_DE_outside_noise")]

# eQTL lane:
#   - non-DE lanes: never eQTL-tested (lane = "not_evaluated")
#   - DE genes:
#       n_cis_total == 0       -> "no_GTEx_data"
#       n_supp_with_repro >= 1 -> "explained"
#       otherwise              -> "unexplained"
m[, eqtl_lane := fcase(
  !(sig_lane %in% c("DE_low", "DE_high")), "not_evaluated",
  n_cis_total == 0L,                       "no_GTEx_data",
  n_supp_with_repro >= 1L,                 "explained",
  default                                = "unexplained")]

# =============================================================================
# STEP 4: Ordered output table
# =============================================================================

cat("\nStep 4: Writing output tables...\n")

setcolorder(m, c(
  "EnsemblID", "ensembl_stable", "Gene_name", "Gene_type", "Chr",
  "baseMean", "raw_log2FC", "raw_FC", "raw_padj",
  "norm_log2FC", "norm_padj",
  "deviation_magnitude", "deviation_vs_cohort_sd",
  "passes_magnitude_filter",
  "low_expr", "high_repeat",
  "sig_lane", "eqtl_lane",
  "n_cis_total", "n_dir_match", "n_supp_with_repro",
  "strongest_supp_pval", "strongest_supp_variant",
  "strongest_dir_pval", "strongest_dir_variant",
  "strongest_overall_pval", "strongest_overall_variant"
))

setorder(m, sig_lane, eqtl_lane, -deviation_magnitude)

fwrite(m, "results/tables/chr21_lane_assignments.csv")
cat("  Wrote results/tables/chr21_lane_assignments.csv\n")

# Lane counts (overall, and after baseMean filter for the headline numbers)
summary_all <- m[, .(n_genes = .N),
                 by = .(sig_lane, eqtl_lane)]
setorder(summary_all, sig_lane, eqtl_lane)
summary_filt <- m[low_expr == FALSE & high_repeat == FALSE,
                  .(n_genes = .N),
                  by = .(sig_lane, eqtl_lane)]
setorder(summary_filt, sig_lane, eqtl_lane)
summary_all[,  scope := "all_chr21"]
summary_filt[, scope := "after_paper_filters"]

lane_summary <- rbindlist(list(summary_all, summary_filt))
fwrite(lane_summary, "results/tables/chr21_lane_summary.csv")
cat("  Wrote results/tables/chr21_lane_summary.csv\n")

# =============================================================================
# STEP 5: Verification
# =============================================================================

cat("\n=== Verification ===\n")
cat(sprintf("Total chr21 genes assigned: %d\n", nrow(m)))
cat(sprintf("Low-expression flagged (q%d):  %d\n",
            as.integer(100 * LOW_EXPR_QUANT), sum(m$low_expr)))
cat(sprintf("High-repeat flagged:          %d\n", sum(m$high_repeat)))
cat(sprintf("Genes with any GTEx eQTL data: %d / %d\n",
            sum(m$n_cis_total > 0), nrow(m)))

cat("\nLane counts (after baseMean q20 + repeat filter):\n")
print(summary_filt)

cat("\nLane counts (all chr21):\n")
print(summary_all)

# Spot checks on canonical genes
cat("\nSpot check (APP, COL18A1, OLIG2, BACE2, MX1, CSTB):\n")
print(m[Gene_name %in% c("APP", "COL18A1", "OLIG2", "BACE2", "MX1", "CSTB"),
        .(Gene_name, baseMean = round(baseMean),
          raw_FC = round(raw_FC, 2),
          norm_log2FC = round(norm_log2FC, 2),
          norm_padj = signif(norm_padj, 2),
          sig_lane, eqtl_lane,
          n_cis_total, n_dir_match, n_supp_with_repro)])

writeLines(capture.output(sessionInfo()),
           "results/tables/chr21_lane_assignment_session_info.txt")

cat("\n=== Lane assignment complete ===\n")
