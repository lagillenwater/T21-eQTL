# 04_chr21_lane_assignment.R
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
#          eQTL lane: cis_eqtl         (gene-level permutation test detects a
#                                       cis-eQTL: q_gene_bh < FDR_GENE - see
#                                       Step 5 in scripts/03_t21_dosage_boxplots.R)
#                     no_cis_eqtl      (permutation test run, did not detect one)
#                     no_GTEx_data     (gene not in GTEx whole-blood
#                                       signif_pairs at all)
#
#          NOTE: eqtl_lane is NOT an "explained by eQTL" claim - see the
#          tight plan (docs/superpowers/plans/2026-08-31-tight-plan.md),
#          which retires that framing, and docs/REPO_STATE.md's decision log.
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
source("scripts/lib/run.R"); run <- load_run(); th <- run$thresholds

cat(sprintf("=== T21-eQTL: Per-Gene Lane Assignment (DE x eQTL) [run: %s] ===\n\n", run$name))

# =============================================================================
# Constants - matched to paper conventions
# =============================================================================

# Thresholds come from the run config (config/runs/<name>.R); local names kept.
ALPHA               <- th$alpha_de          # paper: padj < .01 (Fig. 2B, 3B)
ALPHA_REPRO         <- th$alpha_repro       # within-T21 nominal p, context column only
OUTLIER_FDR         <- th$outlier_fdr       # chr21-internal outlier annotation only
DEVIATION_LFC       <- th$deviation_lfc     # tier 1 (primary, Hunter)
DEVIATION_LFC_T2    <- th$deviation_lfc_t2  # tier 2 (secondary); the cut assign_sig_lane splits on
LOW_EXPR_BASEMEAN   <- th$low_expr_basemean # Hunter et al.'s minimum read coverage
source("scripts/lib/biotypes.R")   # TARGET_BIOTYPES: the chr21 target set, the
                                   # biotypes GTEx whole blood tests

KNOWN_REPEAT_GENES <- c("RPS6KB1", "RPS27", "RPS27L", "RPS27P",
                        "IFNAR1", "IFNAR2", "TPTE", "BAGE", "DAB1")

# =============================================================================
# STEP 1: Load inputs
# =============================================================================

cat("Step 1: Loading inputs...\n")

deseq <- fread(run$table("deseq2_chr21_genes_both_analyses.csv"))
deseq[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
n_before <- nrow(deseq)
deseq <- deseq[Gene_type %in% TARGET_BIOTYPES]
cat(sprintf("  chr21 DESeq2 rows (%s): %d (was %d)\n",
            paste(TARGET_BIOTYPES, collapse = " + "), nrow(deseq), n_before))
print(deseq[, .N, by = Gene_type])

per_var <- fread(run$table("t21_dosage_per_variant.csv"))
cat(sprintf("  Per-variant within-T21 fits: %d\n", nrow(per_var)))

targets <- fread(run$processed("eqtl_target_variants.csv"))
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
# dir_match is ALT-referenced: REF and ALT come from the assembly, not from
# the population, so its meaning flips with an arbitrary label. dir_match_minor
# asks the same question of the MINOR allele (script 03, scripts/lib/alleles.R)
# and is the one comparable across variants and with published eQTL
# directions. Both are locus-level context; neither gates a lane.
per_var[, dir_match_minor := !is.na(gtex_slope_minor) & !is.na(norm_log2FC_obs) &
                             sign(gtex_slope_minor) == sign(norm_log2FC_obs)]
per_var[, supportive_with_repro :=
          dir_match &
          !is.na(t21_slope) & !is.na(t21_p) &
          sign(t21_slope) == sign(gtex_slope) &
          t21_p < ALPHA_REPRO]

# Per-gene aggregation
locus <- per_var[, .(
  n_cis_total       = uniqueN(variant_id),
  n_dir_match       = uniqueN(variant_id[dir_match == TRUE]),
  n_dir_match_minor = uniqueN(variant_id[dir_match_minor == TRUE]),
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
m[is.na(n_dir_match_minor), n_dir_match_minor := 0L]
m[is.na(n_supp_with_repro), n_supp_with_repro := 0L]
m[, raw_FC := 2^raw_log2FC]

# Paper's filters as flags (categorize, do not remove). The low-expression
# cutoff is Hunter et al.'s absolute minimum (matches script 02).
basemean_threshold <- LOW_EXPR_BASEMEAN
cat(sprintf("  baseMean cutoff (absolute): %.2f\n", basemean_threshold))
m[, low_expr   := !is.na(baseMean) & baseMean < basemean_threshold]
m[, high_repeat := Gene_name %in% high_repeat_genes]

m[, deviation_magnitude := abs(norm_log2FC)]

# Significance lane - padj + tier-2 magnitude cut, gated on eligibility.
# The rule (and the reason its fcase order is what it is) lives in
# scripts/lib/lane_rules.R so it can be unit-tested; this call is the only
# place it is applied. It adds eligible_idx, passes_magnitude_filter, sig_lane.
source("scripts/lib/lane_rules.R")
assign_sig_lane(m, alpha = ALPHA, deviation_lfc = DEVIATION_LFC_T2)
# Tier annotation: 1 = Hunter's primary rule, 2 = secondary band, NA = not DE.
m[, tier := fcase(
  sig_lane %in% c("DE_high", "DE_low") & abs(norm_log2FC) >= DEVIATION_LFC, 1L,
  sig_lane %in% c("DE_high", "DE_low"),                                     2L,
  default = NA_integer_)]
cat(sprintf("  Tier 1 (>= %.3f): %d genes; tier 2 (>= %.3f): %d genes\n",
            DEVIATION_LFC, sum(m$tier == 1L, na.rm = TRUE),
            DEVIATION_LFC_T2, sum(m$tier == 2L, na.rm = TRUE)))

cat(sprintf("  Deviating (padj < %.2g AND |corrected log2FC| >= %.3f): %d genes\n",
            ALPHA, DEVIATION_LFC_T2,
            sum(m$sig_lane %in% c("DE_low", "DE_high"))))
cat(sprintf("  Clearing the magnitude cut alone: %d eligible genes\n",
            sum(m$passes_magnitude_filter)))

# Deviation magnitude vs a chr21-internal null. ANNOTATION ONLY - dev_z,
# q_outlier and chr21_k_sensitivity.csv do not gate any lane; the
# assign_sig_lane rule above is the sole classification rule. See docs/REPO_STATE.md decision log.
source("scripts/lib/chr21_threshold.R")

# Null from eligible genes only; z and q reported for all genes so the table
# stays complete. q is NA for ineligible genes - they are caught by the
# High_repeats / Low_expression lanes.
eligible_idx <- m$eligible_idx
null <- chr21_null(m$norm_log2FC[eligible_idx])
cat(sprintf("  chr21 null: center %.4f  MAD %.4f  (n = %d eligible genes)\n",
            null$center, null$scale, null$n))

m[, dev_z := robust_z(norm_log2FC, null)]
m[, q_outlier := NA_real_]
m[eligible_idx, q_outlier := outlier_fdr(dev_z)]

cat(sprintf("  Annotation only - FDR-outlier test at FDR < %.2f flags %d genes (effective k = %.2f)\n",
            OUTLIER_FDR, sum(m$q_outlier < OUTLIER_FDR, na.rm = TRUE),
            effective_k(m$dev_z, m$q_outlier, OUTLIER_FDR)))

sens <- k_sensitivity(m$dev_z[eligible_idx])
fwrite(sens, run$table("chr21_k_sensitivity.csv"))
stopifnot(file.exists(run$table("chr21_k_sensitivity.csv")))
cat("  Wrote", run$table("chr21_k_sensitivity.csv"), "\n")
print(sens)

# eQTL lane, now gated on gene-level permutation significance rather than
# "at least one supportive variant". The old rule scaled with the number of cis
# variants tested (median n_cis 107 for cis_eqtl genes vs 36 for the one
# no_cis_eqtl gene), so it measured variant count more than genetic evidence.
perm <- if (file.exists(run$table("eqtl_gene_level_perm.csv"))) {
  fread(run$table("eqtl_gene_level_perm.csv"))
} else {
  stop("run scripts/03_t21_dosage_boxplots.R --run ", run$name, " first - eqtl_gene_level_perm.csv is missing")
}
m <- merge(m, perm[, .(Gene_name, p_gene_perm, q_gene_bh, cis_eqtl_detected, best_variant)],
           by = "Gene_name", all.x = TRUE)

# Minor-allele reference for the variant the call rests on, so the lane table
# states the direction of the detected eQTL per copy of the allele that is
# rarer in the population rather than per copy of whichever base differs from
# the reference assembly (scripts/lib/alleles.R).
best_allele <- unique(per_var[, .(Gene_name, best_variant = variant_id,
                                  best_minor_allele    = minor_allele,
                                  best_major_allele    = major_allele,
                                  best_maf_gtex        = gtex_maf,
                                  best_maf_gnomad      = gnomad_maf,
                                  best_maf_htp         = htp_maf,
                                  best_minor_concordant = minor_concordant,
                                  best_slope_minor_t21  = t21_slope_minor,
                                  best_slope_minor_gtex = gtex_slope_minor)])
m <- merge(m, best_allele, by = c("Gene_name", "best_variant"), all.x = TRUE)

# eQTL lane:
#   - non-DE lanes: never eQTL-tested (lane = "not_evaluated")
#   - DE genes:
#       n_cis_total == 0                                -> "no_GTEx_data"
#       gene-level permutation q < FDR_GENE (BH)         -> "cis_eqtl"
#       otherwise                                        -> "no_cis_eqtl"
# eqtl_lane is a detection result, not an "explained by eQTL" claim - see
# docs/REPO_STATE.md decision log.
m[, eqtl_lane := fcase(
  !(sig_lane %in% c("DE_low", "DE_high")),             "not_evaluated",
  n_cis_total == 0L,                                   "no_GTEx_data",
  !is.na(cis_eqtl_detected) & cis_eqtl_detected == TRUE, "cis_eqtl",
  default =                                            "no_cis_eqtl")]

# =============================================================================
# STEP 4: Ordered output table
# =============================================================================

cat("\nStep 4: Writing output tables...\n")

setcolorder(m, c(
  "EnsemblID", "ensembl_stable", "Gene_name", "Gene_type", "Chr",
  "baseMean", "raw_log2FC", "raw_FC", "raw_padj",
  "norm_log2FC", "norm_padj",
  "deviation_magnitude", "dev_z", "q_outlier",
  "eligible_idx", "passes_magnitude_filter", "tier",
  "low_expr", "high_repeat",
  "sig_lane", "eqtl_lane",
  "n_cis_total", "n_dir_match", "n_dir_match_minor", "n_supp_with_repro",
  "p_gene_perm", "q_gene_bh", "cis_eqtl_detected",
  "best_minor_allele", "best_major_allele", "best_maf_gtex", "best_maf_gnomad",
  "best_maf_htp", "best_minor_concordant",
  "best_slope_minor_t21", "best_slope_minor_gtex",
  "strongest_supp_pval", "strongest_supp_variant",
  "strongest_dir_pval", "strongest_dir_variant",
  "strongest_overall_pval", "strongest_overall_variant"
))

setorder(m, sig_lane, eqtl_lane, -deviation_magnitude)

fwrite(m, run$table("chr21_lane_assignments.csv"))
cat("  Wrote", run$table("chr21_lane_assignments.csv"), "\n")

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
fwrite(lane_summary, run$table("chr21_lane_summary.csv"))
cat("  Wrote", run$table("chr21_lane_summary.csv"), "\n")

# =============================================================================
# STEP 5: Verification
# =============================================================================

cat("\n=== Verification ===\n")
cat(sprintf("Total chr21 genes assigned: %d\n", nrow(m)))
print(m[, .N, by = Gene_type])
cat(sprintf("Low-expression flagged (baseMean < %g): %d\n",
            LOW_EXPR_BASEMEAN, sum(m$low_expr)))
cat(sprintf("High-repeat flagged:          %d\n", sum(m$high_repeat)))
cat(sprintf("Genes with any GTEx eQTL data: %d / %d\n",
            sum(m$n_cis_total > 0), nrow(m)))

cat("\nLane counts (after baseMean + repeat filter):\n")
print(summary_filt)

cat("\nsig_lane by biotype (all chr21):\n")
print(dcast(m[, .N, by = .(sig_lane, Gene_type)], sig_lane ~ Gene_type,
            value.var = "N", fill = 0L))

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
           run$table("chr21_lane_assignment_session_info.txt"))
stopifnot(file.exists(run$table("chr21_lane_assignment_session_info.txt")))

stopifnot(
  file.exists(run$table("chr21_lane_assignments.csv")),
  file.exists(run$table("chr21_lane_summary.csv")),
  file.exists(run$table("chr21_k_sensitivity.csv"))
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
#
# 2026-08-31  RENAMED eqtl_lane values "explained" -> "cis_eqtl" and
#             "unexplained" -> "no_cis_eqtl"; merged column explained_perm ->
#             cis_eqtl_detected (matches the rename in script 03).
#             "no_GTEx_data" and "not_evaluated" are unchanged. Reason: the
#             tight plan retires any "explained by eQTL" claim - see
#             docs/REPO_STATE.md decision log. Downstream consumers updated:
#             scripts/05_alluvial_lane_assignment.R,
#             scripts/06_chr21_distribution_panel.R,
#             scripts/07_three_panel_figure.R.
#
# 2026-08-31  REPLACED the composition null with a CORRELATION-MATCHED one
#             (partner_null now takes L_ctrl and draws random seed genes, each
#             contributing the median log2FC of its OWN top-20 correlated
#             partners; 300 draws, p reported as (1 + k) / (n_draw + 1)).
#             Reason: the old null drew independent random 20-gene sets, but
#             the observed partners are a co-expression module and move
#             together - measured on this cohort, module medians have SD 0.259
#             against 0.058 for independent sets. Testing a module against an
#             independent null is anti-conservative and drove both PROGRAM-side
#             p-values to exactly 0. Verdicts are unchanged (COL6A1 MIXED
#             p 0.043, OLIG2 MIXED p 0.013, RIPK4 GENE-SPECIFIC p 0.316,
#             TSPEAR GENE-SPECIFIC p 0.492).
#
# 2026-08-31  EXTRACTED the sig_lane rule to scripts/lib/lane_rules.R
#             (assign_sig_lane), with tests/testthat/test-lane-rules.R covering
#             lane reachability, the eligibility-before-magnitude fcase order,
#             the inclusive log2(1.5) boundary, sign routing and NA log2FC.
#             Logic is unchanged - the sig_lane column is byte-identical to the
#             pre-extraction table. Reason: this rule produces the headline
#             classification, lived inline, and had already broken once with an
#             unreachable-lanes regression that only output inspection caught.
#             Side effects: eligible_idx is now a column of the lane table
#             rather than a local vector, and the FDR-outlier print no longer
#             welds the Hunter-rule gene count onto the outlier test's
#             effective k - the two are printed separately, and the outlier
#             line is labelled annotation-only.
# 2026-09-01  ADDED the two-tier scheme: DE lanes now admit tier 2
#             (>= log2(4/3)) as well as tier 1 (>= log2(1.5), Hunter, primary),
#             distinguished by the new `tier` column. passes_magnitude_filter
#             therefore reflects the tier-2 threshold. Composition control and
#             the cis-eQTL permutation run on both tiers.
# 2026-09-04  WIDENED the chr21 target set to protein-coding, lncRNA and
#             pseudogene biotypes (TARGET_BIOTYPES, scripts/lib/biotypes.R):
#             GTEx whole blood tests all of them for cis-eQTLs, so restricting
#             to protein-coding left a third of the eQTL-testable chr21 genes
#             out of the question.
#             REPLACED the q20 baseMean low-expression cutoff with Hunter et
#             al.'s absolute minimum (LOW_EXPR_BASEMEAN = 30), because a
#             quantile of the target set moves when the set changes (25.1 ->
#             4.1 with lncRNA in it). The verification block now prints the
#             sig_lane x biotype table; the lane table already carries
#             Gene_type, so no output schema changed.
# 2026-09-04  RENAMED the composition control to the co-expression
#             neighborhood check (scripts/lib/neighborhood.R) and its labels
#             to the descriptive SHARED / PARTLY-SHARED / NOT-SHARED (were
#             PROGRAM / MIXED / GENE-SPECIFIC). The test is unchanged; the
#             old names claimed a cause (composition, a program) the test
#             does not observe. ADDED partner_z and partner_pctl (the
#             neighborhood shift in SDs and percentile of the matched null)
#             and carried partner_lfc into the lane table. Columns renamed:
#             verdict -> neighborhood, residual_lfc -> gene_minus_partner_lfc,
#             program_share -> share_ratio; output file
#             chr21_composition_control.csv -> chr21_neighborhood_shift.csv.
# 2026-09-10  RENAMED the neighborhood labels SHARED / PARTLY-SHARED /
#             NOT-SHARED to neighbors_shift / neighbors_shift_less /
#             neighbors_no_shift (scripts/lib/neighborhood.R). The rule is
#             unchanged; the new names state the observation (did the
#             partners shift, and by how much relative to the gene) instead
#             of an adjective that needed the legend to decode.
# 2026-09-10  RETIRED the co-expression neighborhood check (Step 3b,
#             scripts/lib/neighborhood.R, chr21_neighborhood_shift.csv, and
#             the lane-table columns neighborhood / partner_lfc / partner_z /
#             gene_minus_partner_lfc), together with scripts 08 and 09.
#             Reason: in whole blood a gene's co-expression partners are
#             largely a cell-type signature, so the check was a proxy for
#             composition; the adjusted run now uses measured CyTOF cell
#             fractions directly, and the neighborhood label had no reading
#             that the cell fractions do not give more plainly. Lane
#             classification and eqtl_lane are unchanged.
