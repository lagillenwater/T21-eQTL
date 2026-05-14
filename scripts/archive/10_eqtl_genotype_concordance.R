# 10_eqtl_genotype_concordance.R
#
# Purpose: Per-variant directional concordance between GTEx whole-blood eQTL
#          slopes and the T21 vs Control allele dosage shift in the HTP cohort.
#          A variant is "concordant" if the dosage shift between T21 and Control
#          (sign of mean(T21) - mean(Control)) times the GTEx slope sign matches
#          the observed sign of the gene's ploidy-normalized log2FC. Aggregates
#          to a per-gene summary.
#
# Inputs:
#   - data/processed/eqtl_target_variants.csv  (GTEx slopes + gene metadata)
#   - data/processed/genotypes_filtered.csv    (long: variant x sample dosage)
#   - data/processed/eqtl_supported_genes.csv  (target gene table)
#
# Outputs:
#   - results/tables/eqtl_genotype_concordance_per_variant.csv
#   - results/tables/eqtl_genotype_concordance_per_gene.csv
#   - results/tables/eqtl_genotype_concordance_session_info.txt
#
# Date: 2026-04-30

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

set.seed(42)

cat("=== T21-eQTL: Genotype-eQTL Directional Concordance ===\n\n")

# =============================================================================
# STEP 1: Load inputs
# =============================================================================

cat("Step 1: Loading filtered genotypes and target variants...\n")

geno    <- fread("data/processed/genotypes_filtered.csv")
targets <- fread("data/processed/eqtl_target_variants.csv")
genes   <- fread("data/processed/eqtl_supported_genes.csv")

stopifnot(all(geno$karyotype %in% c("T21", "Control")))

cat(sprintf("  Genotype rows:         %d\n", nrow(geno)))
cat(sprintf("  Target variant pairs:  %d\n", nrow(targets)))
cat(sprintf("  Target genes:          %d\n", nrow(genes)))

# =============================================================================
# STEP 2: Per-(variant, sample) -> per-variant group means
# =============================================================================

cat("\nStep 2: Computing per-variant T21 vs Control dosage means...\n")

# A subject can sit on a multi-suffix lab_id - collapse to one dosage per
# subject per variant by averaging duplicates (rare; warn if found).
dup_check <- geno[, .N, by = .(variant_id, subject_id)][N > 1]
if (nrow(dup_check) > 0) {
  warning(sprintf("%d (variant,subject) pairs had multiple genotypes; averaged.",
                  nrow(dup_check)))
}

per_subject <- geno[!is.na(alt_dosage),
                    .(alt_dosage = mean(alt_dosage)),
                    by = .(variant_id, subject_id, karyotype)]

# To make T21 (3 alleles) and Control (2 alleles) directly comparable, also
# compute alt allele frequency per subject = alt_dosage / ploidy.
per_subject[, ploidy := fifelse(karyotype == "T21", 3L, 2L)]
per_subject[, alt_freq := alt_dosage / ploidy]

# Per-variant group means and a Welch t-test on alt_freq
per_variant <- per_subject[, .(
    n_T21          = sum(karyotype == "T21"),
    n_Control      = sum(karyotype == "Control"),
    mean_dosage_T21     = mean(alt_dosage[karyotype == "T21"]),
    mean_dosage_Control = mean(alt_dosage[karyotype == "Control"]),
    mean_freq_T21       = mean(alt_freq[karyotype == "T21"]),
    mean_freq_Control   = mean(alt_freq[karyotype == "Control"])
  ),
  by = variant_id
]

# Welch t-test on per-subject alt_freq (vectorized via split-apply)
welch_p <- per_subject[, {
  t21 <- alt_freq[karyotype == "T21"]
  ctl <- alt_freq[karyotype == "Control"]
  if (length(t21) < 2 || length(ctl) < 2 ||
      var(t21) + var(ctl) == 0) {
    .(p_value = NA_real_)
  } else {
    .(p_value = tryCatch(
      t.test(t21, ctl, var.equal = FALSE)$p.value,
      error = function(e) NA_real_))
  }
}, by = variant_id]

per_variant <- merge(per_variant, welch_p, by = "variant_id", all.x = TRUE)
per_variant[, delta_freq := mean_freq_T21 - mean_freq_Control]
per_variant[, delta_dosage := mean_dosage_T21 - mean_dosage_Control]

# =============================================================================
# STEP 3: Join GTEx slope + observed expression direction; score concordance
# =============================================================================

cat("\nStep 3: Scoring directional concordance...\n")

# A target variant can map to multiple genes; keep one row per (variant, gene).
target_slim <- targets[, .(variant_id, ensembl_stable, Gene_name, gene_set,
                           slope, slope_se, pval_nominal, af,
                           raw_log2FC, norm_log2FC, observed_direction)]

scored <- merge(target_slim, per_variant, by = "variant_id",
                allow.cartesian = TRUE)

# Predicted direction of expression deviation from genotype:
#   sign(slope * delta_freq)  (positive => predicts higher expression in T21)
# Concordant if predicted direction == observed_direction (sign of norm_log2FC)
scored[, predicted_direction := sign(slope * delta_freq)]
scored[, concordant := fifelse(
  predicted_direction == 0 | observed_direction == 0,
  NA, predicted_direction == observed_direction)]

cat(sprintf("  Variant-gene rows scored: %d\n", nrow(scored)))
cat(sprintf("  Concordant: %d  Discordant: %d  Indeterminate: %d\n",
            sum(scored$concordant == TRUE,  na.rm = TRUE),
            sum(scored$concordant == FALSE, na.rm = TRUE),
            sum(is.na(scored$concordant))))

scored[, abs_slope := abs(slope)]
setorder(scored, ensembl_stable, -abs_slope)
scored[, abs_slope := NULL]
fwrite(scored, "results/tables/eqtl_genotype_concordance_per_variant.csv")

# =============================================================================
# STEP 4: Per-gene summary
# =============================================================================

cat("\nStep 4: Aggregating to per-gene summary...\n")

per_gene <- scored[, .(
    gene_set      = first(gene_set),
    raw_log2FC    = first(raw_log2FC),
    norm_log2FC   = first(norm_log2FC),
    observed_direction = first(observed_direction),
    n_variants    = .N,
    n_concordant  = sum(concordant == TRUE,  na.rm = TRUE),
    n_discordant  = sum(concordant == FALSE, na.rm = TRUE),
    n_indet       = sum(is.na(concordant)),
    frac_concordant = mean(concordant, na.rm = TRUE),
    mean_slope    = mean(slope),
    mean_delta_freq = mean(delta_freq, na.rm = TRUE),
    median_p      = median(p_value, na.rm = TRUE),
    n_sig_dosage  = sum(p_value < 0.05, na.rm = TRUE)
  ),
  by = .(ensembl_stable, Gene_name)
]

per_gene[, gene_supported := frac_concordant >= 0.5 & n_concordant >= 1]

setorder(per_gene, gene_set, -frac_concordant, -n_variants)
fwrite(per_gene, "results/tables/eqtl_genotype_concordance_per_gene.csv")

# =============================================================================
# STEP 5: Verification summary
# =============================================================================

cat("\n=== Summary ===\n")
print(per_gene[, .(
  n_genes            = .N,
  n_supported       = sum(gene_supported, na.rm = TRUE),
  n_not_supported   = sum(!gene_supported, na.rm = TRUE),
  median_frac_concord = median(frac_concordant, na.rm = TRUE)
), by = gene_set])

writeLines(capture.output(sessionInfo()),
           "results/tables/eqtl_genotype_concordance_session_info.txt")

cat("\n=== Concordance analysis complete ===\n")
