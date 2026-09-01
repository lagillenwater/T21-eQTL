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
DEVIATION_LFC       <- log2(1.5)   # tier 1 (primary, Hunter)
DEVIATION_LFC_T2    <- log2(4/3)   # tier 2 (secondary)   # Hunter et al.'s FC >= 1.5 cut, applied on
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

m[, deviation_magnitude := abs(norm_log2FC)]

# Significance lane - Hunter's padj + 1.5-fold rule, gated on eligibility.
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

cat(sprintf("  Hunter rule (padj < %.2g AND |corrected log2FC| >= %.3f): %d genes\n",
            ALPHA, DEVIATION_LFC,
            sum(m$sig_lane %in% c("DE_low", "DE_high"))))
cat(sprintf("  Clearing the magnitude cut alone: %d eligible genes\n",
            sum(m$passes_magnitude_filter)))

# Deviation magnitude vs a chr21-internal null. ANNOTATION ONLY - dev_z,
# q_outlier and chr21_k_sensitivity.csv do not gate any lane; the Hunter rule
# above is the sole classification rule. See docs/REPO_STATE.md decision log.
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
fwrite(sens, "results/tables/chr21_k_sensitivity.csv")
stopifnot(file.exists("results/tables/chr21_k_sensitivity.csv"))
cat("  Wrote results/tables/chr21_k_sensitivity.csv\n")
print(sens)

# eQTL lane, now gated on gene-level permutation significance rather than
# "at least one supportive variant". The old rule scaled with the number of cis
# variants tested (median n_cis 107 for cis_eqtl genes vs 36 for the one
# no_cis_eqtl gene), so it measured variant count more than genetic evidence.
perm <- if (file.exists("results/tables/eqtl_gene_level_perm.csv")) {
  fread("results/tables/eqtl_gene_level_perm.csv")
} else {
  stop("run scripts/03_t21_dosage_boxplots.R first - eqtl_gene_level_perm.csv is missing")
}
m <- merge(m, perm[, .(Gene_name, p_gene_perm, q_gene_bh, cis_eqtl_detected, best_variant)],
           by = "Gene_name", all.x = TRUE)

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
# STEP 3b: Composition control for deviating (DE_low / DE_high) genes
# =============================================================================
#
# A chr21 gene can appear to deviate from dosage expectation because the
# blood cell type that expresses it changed in abundance in T21, not because
# of any regulatory effect on the gene itself. Composition shifts move whole
# co-expression programs, not single genes: for each deviating gene, find its
# 20 most co-expressed non-chr21 genes using CONTROLS ONLY (so the karyotype
# effect cannot leak into the correlation), take the median T21-vs-Control
# log2FC of those partners, and compare it to a CORRELATION-MATCHED null -
# random seed genes from the same pool, each contributing the median log2FC of
# its own top-20 correlated partners. Matching matters: a co-expression module
# moves together, so its median log2FC is several times more variable than an
# independent 20-gene set's, and an independent null would call almost any
# partner shift significant. If the partners shift with the gene, that looks
# like a program (composition or shared pathway), not gene-specific dosage
# regulation.

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

null <- partner_null(L_ctrl, lfc_bg, n_partners = 20, n_draw = 300, seed = 1)
cat(sprintf("  Correlation-matched null: %d module draws, median %.4f, SD %.4f\n",
            length(null), median(null), sd(null)))

composition <- rbindlist(lapply(deviating_genes, function(g) {
  gene_lfc <- genome_lfc[[g]]
  gene_ctrl <- L[g, karyotype == "Control"]
  partner_lfc <- partner_shift(g, L_ctrl, gene_ctrl, lfc_bg, n_partners = 20)
  p_partners <- partner_p(gene_lfc, partner_lfc, null)
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
  "eligible_idx", "passes_magnitude_filter", "tier",
  "low_expr", "high_repeat",
  "sig_lane", "eqtl_lane",
  "verdict", "residual_lfc",
  "n_cis_total", "n_dir_match", "n_supp_with_repro",
  "p_gene_perm", "q_gene_bh", "cis_eqtl_detected",
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
