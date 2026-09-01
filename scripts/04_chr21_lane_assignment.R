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
#   - results/tables/deseq2_chr21_genes_both_analyses.csv
#   - results/tables/t21_dosage_per_variant.csv
#   - data/processed/eqtl_target_variants.csv
#   - data/processed/blacklisted_genes.csv  (high-repeat flag)
#
# Outputs:
#   - results/tables/chr21_lane_assignments.csv
#   - results/tables/chr21_lane_summary.csv
#   - results/tables/chr21_k_sensitivity.csv
#   - results/tables/chr21_lane_assignment_session_info.txt
#
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
OUTLIER_FDR         <- 0.10    # FDR-controlled robust outlier test against a
                               # chr21-internal median/MAD null. Genes that do
                               # not clear this are flagged "Expected_dosage"
                               # - their statistically significant deviations
                               # are within typical chr21 variation, so no
                               # eQTL explanation is sought. This filter is
                               # also applied upstream in script 02 (gene
                               # selection) so failing genes have no genotype
                               # data. Eligibility (low_expr / high_repeat)
                               # also determines lane routing directly: those
                               # genes are excluded from the null because they
                               # are too noisy to assess, not because they
                               # follow expected dosage, so the sig_lane
                               # fcase below must test high_repeat/low_expr
                               # BEFORE passes_magnitude_filter.
DEVIATION_LFC       <- log2(1.5)   # Hunter et al.'s FC >= 1.5 cut, applied on
                                   # the ploidy-corrected log2FC scale
LOW_EXPR_QUANT      <- 0.20    # paper: "second quintile of baseMean"
RESTRICT_TO_PROTEIN_CODING <- TRUE   # restrict chr21 set + cohort-noise
                                     # reference to protein-coding genes

KNOWN_REPEAT_GENES <- c("RPS6KB1", "RPS27", "RPS27L", "RPS27P",
                        "IFNAR1", "IFNAR2", "TPTE", "BAGE", "DAB1")

# =============================================================================
# STEP 1: Load inputs
# =============================================================================

cat("Step 1: Loading inputs...\n")

deseq <- fread("results/tables/deseq2_chr21_genes_both_analyses.csv")
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

# Deviation magnitude vs a chr21-internal null (computed BEFORE lane
# assignment so the outlier filter can gate the eQTL lane downstream).
source("scripts/lib/chr21_threshold.R")

m[, deviation_magnitude := abs(norm_log2FC)]

# Null from eligible genes only; z and q reported for all genes so the table
# stays complete. q is NA for ineligible genes, which fails the filter by
# construction - they are caught by the High_repeats / Low_expression lanes.
eligible_idx <- !m$low_expr & !m$high_repeat & !is.na(m$norm_log2FC)
null <- chr21_null(m$norm_log2FC[eligible_idx])
cat(sprintf("  chr21 null: center %.4f  MAD %.4f  (n = %d eligible genes)\n",
            null$center, null$scale, null$n))

m[, dev_z := robust_z(norm_log2FC, null)]
m[, q_outlier := NA_real_]
m[eligible_idx, q_outlier := outlier_fdr(dev_z)]
m[, passes_magnitude_filter := eligible_idx & abs(norm_log2FC) >= DEVIATION_LFC]

cat(sprintf("  Outlier test at FDR < %.2f: %d genes (effective k = %.2f)\n",
            OUTLIER_FDR, sum(m$passes_magnitude_filter),
            effective_k(m$dev_z, m$q_outlier, OUTLIER_FDR)))

sens <- k_sensitivity(m$dev_z[eligible_idx])
fwrite(sens, "results/tables/chr21_k_sensitivity.csv")
stopifnot(file.exists("results/tables/chr21_k_sensitivity.csv"))
cat("  Wrote results/tables/chr21_k_sensitivity.csv\n")
print(sens)

# Significance lane - eligibility (high_repeat / low_expr) is tested FIRST,
# ahead of the outlier filter: those genes were excluded from the null
# because they are too noisy to assess, not because they follow expected
# dosage, so labeling them Expected_dosage would be a false claim. Only
# eligible genes fall through to the outlier filter; genes whose
# ploidy-corrected deviation is within typical chr21 variation are
# categorized as Expected_dosage regardless of raw FC sign or padj
# (statistical significance at large n is not the same as biological
# departure from ploidy expectation). Surviving genes get the standard
# categorization.
m[, sig_lane := fcase(
  high_repeat == TRUE,                                     "High_repeats",
  low_expr == TRUE,                                        "Low_expression",
  passes_magnitude_filter == FALSE,                        "Expected_dosage",
  !is.na(norm_padj) & norm_padj < ALPHA & norm_log2FC < 0, "DE_low",
  !is.na(norm_padj) & norm_padj < ALPHA & norm_log2FC > 0, "DE_high",
  default                                                = "Not_DE_outside_noise")]

# eQTL lane, now gated on gene-level permutation significance rather than
# "at least one supportive variant". The old rule scaled with the number of cis
# variants tested (median n_cis 107 for explained genes vs 36 for the one
# unexplained gene), so it measured variant count more than genetic evidence.
perm <- if (file.exists("results/tables/eqtl_gene_level_perm.csv")) {
  fread("results/tables/eqtl_gene_level_perm.csv")
} else {
  stop("run scripts/03_t21_dosage_boxplots.R first - eqtl_gene_level_perm.csv is missing")
}
m <- merge(m, perm[, .(Gene_name, p_gene_perm, q_gene_bh, explained_perm, best_variant)],
           by = "Gene_name", all.x = TRUE)

# eQTL lane:
#   - non-DE lanes: never eQTL-tested (lane = "not_evaluated")
#   - DE genes:
#       n_cis_total == 0                                -> "no_GTEx_data"
#       gene-level permutation q < FDR_GENE (BH)         -> "explained"
#       otherwise                                        -> "unexplained"
m[, eqtl_lane := fcase(
  !(sig_lane %in% c("DE_low", "DE_high")),         "not_evaluated",
  n_cis_total == 0L,                               "no_GTEx_data",
  !is.na(explained_perm) & explained_perm == TRUE, "explained",
  default =                                        "unexplained")]

# =============================================================================
# STEP 3b: Composition control for deviating (DE_low / DE_high) genes
# =============================================================================
#
# A chr21 gene can appear to deviate from dosage expectation because the
# blood cell type that expresses it changed in abundance in T21, not because
# of any regulatory effect on the gene itself. Composition shifts move whole
# co-expression programs, not single genes: for each deviating gene, find its
# 20 most co-expressed non-chr21 genes using CONTROLS ONLY (so the karyotype
# effect cannot leak into the correlation), take the median T21-vs-Control
# log2FC of those partners, and compare it to a null of 2000 random 20-gene
# sets. If the partners shift with the gene, that looks like a program
# (composition or shared pathway), not gene-specific dosage regulation.

cat("\nStep 3b: Composition control for deviating genes...\n")

source("scripts/lib/composition.R")

count_mat <- fread("data/processed/count_matrix.csv")
cohort <- fread("data/processed/analysis_cohort.csv")
de_all <- fread("results/tables/deseq2_all_genes_ploidy_normalized.csv")

lab_ids <- intersect(setdiff(names(count_mat), c("EnsemblID", "Gene_name", "Chr")),
                     cohort$LabID)
karyotype <- cohort$Karyotype[match(lab_ids, cohort$LabID)]

M <- as.matrix(count_mat[, ..lab_ids])
rownames(M) <- count_mat$Gene_name
L <- log2(t(t(M) / colSums(M) * 1e6) + 1)

partner_pool <- rowMeans(M) >= 25 & count_mat$Chr != "chr21"
L_ctrl <- L[partner_pool, karyotype == "Control"]

genome_lfc <- setNames(de_all$log2FoldChange, de_all$Gene_name)
lfc_bg <- genome_lfc[rownames(L_ctrl)]
lfc_bg <- lfc_bg[!is.na(lfc_bg)]

deviating_genes <- m$Gene_name[m$sig_lane %in% c("DE_high", "DE_low")]
deviating_genes <- intersect(deviating_genes, rownames(L))

null <- partner_null(lfc_bg, n_partners = 20, n_draw = 2000, seed = 1)

composition <- rbindlist(lapply(deviating_genes, function(g) {
  gene_lfc <- genome_lfc[[g]]
  gene_ctrl <- L[g, karyotype == "Control"]
  partner_lfc <- partner_shift(g, L_ctrl, gene_ctrl, lfc_bg, n_partners = 20)
  s <- sign(gene_lfc)
  p_partners <- mean(s * null >= s * partner_lfc)
  verdict <- composition_verdict(gene_lfc, partner_lfc, p_partners)
  data.table(
    Gene_name    = g,
    gene_lfc     = gene_lfc,
    partner_lfc  = partner_lfc,
    p_partners   = p_partners,
    program_share = partner_lfc / gene_lfc,
    residual_lfc = gene_lfc - partner_lfc,
    verdict      = verdict
  )
}))

fwrite(composition, "results/tables/chr21_composition_control.csv")
stopifnot(file.exists("results/tables/chr21_composition_control.csv"))
cat("  Wrote results/tables/chr21_composition_control.csv\n")
print(composition)

m <- merge(m, composition[, .(Gene_name, verdict, residual_lfc)],
           by = "Gene_name", all.x = TRUE)

# =============================================================================
# STEP 4: Ordered output table
# =============================================================================

cat("\nStep 4: Writing output tables...\n")

setcolorder(m, c(
  "EnsemblID", "ensembl_stable", "Gene_name", "Gene_type", "Chr",
  "baseMean", "raw_log2FC", "raw_FC", "raw_padj",
  "norm_log2FC", "norm_padj",
  "deviation_magnitude", "dev_z", "q_outlier",
  "passes_magnitude_filter",
  "low_expr", "high_repeat",
  "sig_lane", "eqtl_lane",
  "verdict", "residual_lfc",
  "n_cis_total", "n_dir_match", "n_supp_with_repro",
  "p_gene_perm", "q_gene_bh", "explained_perm",
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

stopifnot(
  file.exists("results/tables/chr21_lane_assignments.csv"),
  file.exists("results/tables/chr21_lane_summary.csv"),
  file.exists("results/tables/chr21_k_sensitivity.csv")
)

cat("\n=== Lane assignment complete ===\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  REPLACED the cohort-SD magnitude filter with the chr21-internal
#             FDR outlier test; dropped column deviation_vs_cohort_sd in favour
#             of dev_z and q_outlier; added chr21_k_sensitivity.csv.
#             Reason: ploidy normalization does not act on diploid genes
#             (mean |raw - norm| 0.0048 off chr21 vs 0.583 on it), so their
#             spread measured a different quantity; and a 1-SD cut selects the
#             top ~third of any distribution (18.3% of non-chr21 genes cleared
#             it themselves, vs 20.6% of chr21 - binomial p = 0.26). The null is
#             now estimated AFTER the expression and repeat filters, because
#             log2FC variance scales with counts.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
#
# 2026-08-31  REPLACED the eqtl_lane rule "at least one cis variant is
#             direction-matched and reproduces at t21_p < 0.05" with gene-level
#             permutation significance at BH FDR < 0.05.
#             Reason: the old rule scaled with the number of cis variants
#             tested and had no multiplicity control. The current three-gene
#             DE set spans 21 to 136 cis variants per gene (TSPEAR 21, OLIG2
#             90, COL6A1 136). An earlier, larger pre-correction DE gene set
#             (no longer current) spanned 21 to 1083 cis variants per gene,
#             with median n_cis 107 for explained genes vs 36 for the one
#             unexplained gene, and RBM11 called explained on 1 supporting
#             variant of 83 where chance predicts ~4 - the clearest evidence
#             the old rule was measuring variant count, not genetic evidence.
#             The lead variant was rejected as an alternative because in LD it
#             is frequently a tag, not the causal variant.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
#
# 2026-08-31  REPLACED the FDR-outlier test driving passes_magnitude_filter
#             with Hunter et al.'s own classification: deviating = padj < 0.01
#             AND abs(norm_log2FC) >= log2(1.5) (DEVIATION_LFC), gated on
#             eligibility (eligible_idx) so it is never NA. dev_z / q_outlier
#             are retained as annotation columns; the sig_lane fcase order is
#             unchanged.
#             Spec: docs/superpowers/plans/2026-08-31-tight-plan.md (Task A)
#
# 2026-08-31  ADDED a composition control for DE_high/DE_low genes (Step 3b):
#             for each deviating gene, find its 20 most co-expressed non-chr21
#             genes using CONTROLS ONLY, take their median T21-vs-Control
#             log2FC, and compare it to a null of 2000 random 20-gene sets.
#             Reason: a chr21 gene can appear to deviate because the blood
#             cell type expressing it changed abundance in T21, not because
#             of a regulatory effect on the gene - composition shifts move
#             whole co-expression programs, not single genes. Writes
#             results/tables/chr21_composition_control.csv and merges
#             verdict / residual_lfc into chr21_lane_assignments.csv.
#             New library: scripts/lib/composition.R.
#             Spec: docs/superpowers/plans/2026-08-31-tight-plan.md (Task B)
