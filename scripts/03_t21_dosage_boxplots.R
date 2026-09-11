# 11_t21_dosage_boxplots.R
#
# Purpose: Within-T21 cis-eQTL test. For each target gene, regress expression
#          on alt-allele dosage in T21 subjects only (per variant), run the
#          gene-level permutation test that classifies the gene, and run the
#          standalone negative and positive controls of that test. The
#          representative supportive variant per gene is kept as a legacy
#          table.
#
# Inputs:
#   - data/processed/genotypes_filtered.csv
#   - data/processed/eqtl_target_variants.csv
#   - data/processed/eqtl_supported_genes.csv
#   - data/processed/count_matrix.csv
#   - data/processed/sample_metadata.csv
#
# Outputs:
#   - results/tables/t21_dosage_per_variant.csv       (deviating genes only)
#   - results/tables/eqtl_gene_level_perm.csv          (the classification test)
#   - results/tables/eqtl_control_negative.csv         (unlinked-variant decoys)
#   - results/tables/eqtl_control_positive.csv         (strong GTEx eGenes)
#   - results/tables/eqtl_controls_summary.csv
#   - results/tables/t21_representative_variants.csv (legacy, archive/14)
#   - results/tables/t21_dosage_session_info.txt
#   Per-gene dosage panels are drawn by scripts/11_eqtl_figures.R from these
#   tables, so this script writes no figure.
#
# Date: 2026-04-30

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source("scripts/lib/eqtl_fit.R")
source("scripts/lib/eqtl_controls.R")
source("scripts/lib/run.R"); run <- load_run()

set.seed(42)

ALPHA_REPRO <- run$thresholds$alpha_repro   # within-T21 nominal p for a cis variant to
                                            # count as reproducible (context column only)

cat(sprintf("=== T21-eQTL: Within-T21 Dosage-Expression Boxplots [run: %s] ===\n\n", run$name))

# =============================================================================
# STEP 1: Load inputs and restrict to T21
# =============================================================================

cat("Step 1: Loading inputs and restricting to T21...\n")

geno    <- fread(run$processed("genotypes_filtered.csv"))
targets <- fread(run$processed("eqtl_target_variants.csv"))
genes   <- fread(run$processed("eqtl_supported_genes.csv"))
meta    <- run$cohort()
# Expression: the run's artifact (log2-CPM, covariate effects removed with
# karyotype kept; script 01). Baseline: plain log2-CPM.
E       <- run$expression()

geno_t21 <- geno[karyotype == "T21" & !is.na(alt_dosage)]
meta_t21 <- meta[Karyotype == "T21"]
meta_t21[, subject_id := sub("[A-Z][0-9]*$", "", LabID)]

# Subject -> LabID used in count matrix (T21 only)
subj_to_labid <- setNames(meta_t21$LabID, meta_t21$subject_id)

cat(sprintf("  T21 subjects with genotypes: %d\n",
            uniqueN(geno_t21$subject_id)))
cat(sprintf("  T21 subjects in metadata:    %d\n",
            uniqueN(meta_t21$subject_id)))
shared_subjects <- intersect(geno_t21$subject_id, meta_t21$subject_id)
cat(sprintf("  Shared (will be analyzed):   %d\n", length(shared_subjects)))

# Collapse subject duplicates by averaging dosage (rare)
geno_t21 <- geno_t21[subject_id %in% shared_subjects,
                     .(alt_dosage = mean(alt_dosage)),
                     by = .(variant_id, subject_id)]

# =============================================================================
# STEP 2: Long-format expression for the 82 target genes (T21 only)
# =============================================================================

cat("\nStep 2: Building long-format expression matrix...\n")

target_names <- unique(targets$Gene_name)
t21_labids   <- intersect(colnames(E), subj_to_labid[shared_subjects])
E_t <- E[rownames(E) %in% target_names, t21_labids, drop = FALSE]
# Rows are keyed by gene name; a target symbol carried by two annotation rows
# would be silently resolved to the first, so refuse that case outright.
dup_targets <- unique(rownames(E_t)[duplicated(rownames(E_t))])
if (length(dup_targets)) stop("target gene name(s) duplicated in the expression artifact: ", paste(dup_targets, collapse = ", "))
cat(sprintf("  Target genes matched in the expression artifact: %d / %d\n",
            nrow(E_t), length(target_names)))

labid_to_subj <- setNames(meta_t21$subject_id, meta_t21$LabID)
expr_long <- data.table(Gene_name = rep(rownames(E_t), times = ncol(E_t)),
                        LabID     = rep(colnames(E_t), each = nrow(E_t)),
                        expr      = as.vector(E_t))
expr_long[, subject_id := labid_to_subj[LabID]]
expr_long <- merge(expr_long, unique(targets[, .(Gene_name, ensembl_stable)]),
                   by = "Gene_name")

cat(sprintf("  Genes with expression rows:   %d\n", uniqueN(expr_long$Gene_name)))
cat(sprintf("  Subjects in expression long:  %d\n", uniqueN(expr_long$subject_id)))

# =============================================================================
# STEP 3: Per-(variant, gene) within-T21 regression of expression on dosage
# =============================================================================

cat("\nStep 3: Regressing expression ~ dosage in T21 (per variant, per gene)\n")

# variant -> gene mapping (a variant can map to multiple target genes)
var_gene <- unique(targets[, .(variant_id, ensembl_stable, Gene_name,
                               gene_set, slope, slope_se, pval_nominal,
                               raw_log2FC, norm_log2FC, observed_direction)])
setnames(var_gene, "slope",        "gtex_slope")
setnames(var_gene, "slope_se",     "gtex_slope_se")
setnames(var_gene, "pval_nominal", "gtex_pval")

# Join genotypes to expression by subject, then to gene per variant mapping
geno_expr <- merge(
  geno_t21,
  expr_long[, .(ensembl_stable, subject_id, expr)],
  by = "subject_id", allow.cartesian = TRUE
)

geno_expr <- merge(geno_expr, var_gene,
                   by = c("variant_id", "ensembl_stable"),
                   all = FALSE)

cat(sprintf("  (variant, gene, subject) rows: %d\n", nrow(geno_expr)))

# Per-(variant, gene) linear regression via scripts/lib/eqtl_fit.R's
# vectorized closed-form fit (single-column call here; Tasks 6/7 reuse the
# same function across all variants of a gene at once). The n < 5 guard is
# preserved explicitly: fit_variants() itself only errors below n = 3, so
# the stricter minimum-sample-size behaviour of this script would silently
# loosen without this check.
fit_table <- geno_expr[, {
  x <- alt_dosage
  y <- expr
  n <- .N
  if (n < 5 || var(x) == 0) {
    .(t21_n = n, t21_slope = NA_real_, t21_se = NA_real_,
      t21_t = NA_real_, t21_p = NA_real_)
  } else {
    fit <- fit_variants(matrix(x, ncol = 1), y)
    .(t21_n     = n,
      t21_slope = fit$slope[1],
      t21_se    = fit$se[1],
      t21_t     = fit$t[1],
      t21_p     = fit$p[1])
  }
}, by = .(variant_id, ensembl_stable, Gene_name, gene_set,
          gtex_slope, gtex_pval, observed_direction)]

fit_table[, supportive := !is.na(t21_slope) &
                          sign(t21_slope) == sign(gtex_slope) &
                          sign(t21_slope) != 0]

# Positive-control genes (script 02) ride along in the fits for the control
# test below but stay out of every main-result table: script 04 reads
# t21_dosage_per_variant.csv for its locus-level columns.
fit_pos   <- fit_table[gene_set == "positive_control"]
fit_table <- fit_table[gene_set != "positive_control"]

cat(sprintf("  Variants tested: %d  Supportive: %d\n",
            nrow(fit_table), sum(fit_table$supportive, na.rm = TRUE)))

fwrite(fit_table, run$table("t21_dosage_per_variant.csv"))
stopifnot(file.exists(run$table("t21_dosage_per_variant.csv")))

# =============================================================================
# GENE-LEVEL PERMUTATION SIGNIFICANCE (GTEx / FastQTL eGene procedure)
# =============================================================================
# The "any of N variants" rule mostly measures N (21 to 136 here). Permuting
# expression labels and taking the smallest p across variants per permutation
# builds the null of the BEST variant, handling multiplicity and LD together.
# Unlike picking the lead variant this does not assume the top association is
# causal - in LD the lead variant is frequently only a tag.
# The runner (scripts/lib/eqtl_controls.R::gene_level_tests) is shared with
# the two standalone controls below, so all three sets get the identical test.

N_PERM   <- as.integer(Sys.getenv("T21_N_PERM", "1000"))
FDR_GENE <- run$thresholds$fdr_gene
DECOY_MIN_DISTANCE <- run$thresholds$decoy_min_distance
cat(sprintf("\n=== Gene-level permutation (%d permutations) ===\n", N_PERM))

gwide <- dcast(geno_t21, subject_id ~ variant_id, value.var = "alt_dosage")
subj  <- gwide$subject_id
G_all <- as.matrix(gwide[, -1])
de_genes    <- unique(fit_table$Gene_name)
variants_of <- split(fit_table$variant_id, fit_table$Gene_name)
expr_of     <- function(g) expr_from_matrix(g, subj, E, meta_t21)

perm_res <- gene_level_tests(de_genes, variants_of, G_all, expr_of,
                             n_perm = N_PERM, seed_base = 2026L, fdr = FDR_GENE)
setnames(perm_res, "detected", "cis_eqtl_detected")
setorder(perm_res, p_gene_perm)
print(as.data.frame(perm_res))

fwrite(perm_res, run$table("eqtl_gene_level_perm.csv"))
if (!file.exists(run$table("eqtl_gene_level_perm.csv"))) {
  stop("failed to write ", run$table("eqtl_gene_level_perm.csv"))
}
cat("  Wrote", run$table("eqtl_gene_level_perm.csv"), "\n")

# =============================================================================
# STANDALONE CONTROLS for the gene-level test
# =============================================================================
# Negative: each deviating gene's expression against the cis variant set of
# another tested gene at least DECOY_MIN_DISTANCE away (no LD with its own
# locus), choosing the candidate with the closest variant count so the decoy
# test carries the same multiplicity. Real genotypes, the same subjects, the
# same test; detections should sit near the FDR level.
# Positive: the strongest GTEx whole-blood eGenes among expressed, non-repeat,
# non-deviating chr21 genes (script 02, gene_set == "positive_control"),
# tested on their own cis variants. Most should be detected; if not, the test
# lacks power here and "tested, not detected" carries no weight.

cat("\n=== Standalone controls ===\n")

tss_vec  <- genes[Gene_name %in% de_genes & !is.na(tss), setNames(tss, Gene_name)]
nvar_vec <- vapply(variants_of[names(tss_vec)], length, integer(1))
decoys   <- assign_decoys(tss_vec, nvar_vec, DECOY_MIN_DISTANCE)
neg_genes      <- decoys[!is.na(decoy_gene), Gene_name]
decoy_variants <- setNames(lapply(neg_genes, function(g)
  variants_of[[decoys[Gene_name == g, decoy_gene]]]), neg_genes)
neg_res <- gene_level_tests(neg_genes, decoy_variants, G_all, expr_of,
                            n_perm = N_PERM, seed_base = 3026L, fdr = FDR_GENE)
neg_res <- merge(decoys, neg_res, by = "Gene_name", all.x = TRUE)
neg_res[is.na(detected), detected := FALSE]
setorder(neg_res, p_gene_perm, na.last = TRUE)
fwrite(neg_res, run$table("eqtl_control_negative.csv"))
cat(sprintf("  Negative (unlinked decoy variants): %d of %d detected at q < %.2f\n",
            sum(neg_res$detected), sum(!is.na(neg_res$p_gene_perm)), FDR_GENE))
print(as.data.frame(neg_res[, .(Gene_name, decoy_gene, distance_mb = round(distance / 1e6, 1),
                                n_variants, p_gene_perm, q_gene_bh, detected)]))

pos_genes    <- genes[gene_set == "positive_control", Gene_name]
pos_variants <- split(fit_pos$variant_id, fit_pos$Gene_name)
pos_res <- gene_level_tests(pos_genes, pos_variants, G_all, expr_of,
                            n_perm = N_PERM, seed_base = 4026L, fdr = FDR_GENE)
pos_res <- merge(genes[gene_set == "positive_control",
                       .(Gene_name, gtex_min_p, norm_log2FC, norm_padj)],
                 pos_res, by = "Gene_name")
setorder(pos_res, gtex_min_p)
fwrite(pos_res, run$table("eqtl_control_positive.csv"))
cat(sprintf("  Positive (strong GTEx eGenes): %d of %d detected at q < %.2f\n",
            sum(pos_res$detected), sum(!is.na(pos_res$p_gene_perm)), FDR_GENE))
print(as.data.frame(pos_res[, .(Gene_name, gtex_min_p = signif(gtex_min_p, 2),
                                n_variants, p_gene_perm, q_gene_bh, detected)]))

ctrl_summary <- data.table(
  set         = c("observed_deviating", "negative_unlinked_variants", "positive_gtex_egenes"),
  n_tested    = c(sum(!is.na(perm_res$p_gene_perm)), sum(!is.na(neg_res$p_gene_perm)),
                  sum(!is.na(pos_res$p_gene_perm))),
  n_detected  = c(sum(perm_res$cis_eqtl_detected), sum(neg_res$detected), sum(pos_res$detected)),
  expectation = c("the result", sprintf("about %.0f%% (the FDR level)", 100 * FDR_GENE),
                  "most detected"))
ctrl_summary[, pct_detected := round(100 * n_detected / n_tested, 1)]
fwrite(ctrl_summary, run$table("eqtl_controls_summary.csv"))
stopifnot(file.exists(run$table("eqtl_control_negative.csv")),
          file.exists(run$table("eqtl_control_positive.csv")),
          file.exists(run$table("eqtl_controls_summary.csv")))
print(as.data.frame(ctrl_summary))

# =============================================================================
# STEP 4: Pick representative variant per gene (most sig supportive eQTL)
# =============================================================================

cat("\nStep 4: Selecting representative variant per gene...\n")

supportive <- fit_table[supportive == TRUE & !is.na(t21_p)]
setorder(supportive, ensembl_stable, t21_p)
representatives <- supportive[, head(.SD, 1L), by = ensembl_stable]

fwrite(representatives,
       run$table("t21_representative_variants.csv"))
stopifnot(file.exists(run$table("t21_representative_variants.csv")))

n_rep_low  <- sum(representatives$gene_set == "DE_low_FC")
n_rep_high <- sum(representatives$gene_set == "Sig_high_FC")
cat(sprintf("  Representative-variant genes:  %d total ",
            nrow(representatives)))
cat(sprintf("(DE_low_FC: %d / 66, Sig_high_FC: %d / 16)\n",
            n_rep_low, n_rep_high))

# =============================================================================
# STEP 6: Verification summary
# =============================================================================

cat("\n=== Summary ===\n")
print(representatives[, .(
  n_genes = .N,
  median_t21_p = median(t21_p, na.rm = TRUE),
  median_n_subj = median(t21_n)
), by = gene_set])

writeLines(capture.output(sessionInfo()),
           run$table("t21_dosage_session_info.txt"))
stopifnot(file.exists(run$table("t21_dosage_session_info.txt")))

cat("\n=== Script 03 complete ===\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  Replaced the inline closed-form regression duplicate in the
#             fit_table block with a call to fit_variants() in
#             scripts/lib/eqtl_fit.R. The inline code was not an lm() call -
#             it was already the same closed-form slope/se/t/p math, just
#             copy-pasted in this script. Consolidated so there is a single
#             tested implementation instead of two copies that could drift,
#             since Tasks 6 and 7 need the same fit_variants/perm_min_p/
#             gene_level_p primitives for gene-level permutation testing.
#             Verified numerically identical to the previous output
#             (max |slope diff| ~ 9.9e-14, max |p diff| ~ 1.0e-15, both well
#             below the 1e-8 tolerance). The n < 5 sample-size guard from the
#             original inline code is preserved explicitly around the
#             fit_variants() call.
#             Reason: single tested regression implementation; Tasks 6/7
#             reuse.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
#
# 2026-08-31  ADDED negative controls C1 (genotype permutation) and C2
#             (direction flip) -> results/tables/eqtl_negative_controls.csv.
#             Reason: the "explained" call had no null. With the current
#             3-gene DE set (OLIG2, COL6A1, TSPEAR; 247 variants), the
#             observed and direction-flipped rates are both 3/3 (100.0%),
#             and the genotype-permuted rate (mean of 20 shuffles) is
#             15.0% (per-shuffle range 0.0-66.7%) -- well below the
#             observed rate but nonzero, consistent with 3 genes giving
#             only a 0/33/67/100% resolution ladder. The "any variant"
#             rule is not fully discriminating at this gene count: a
#             larger DE gene set would be needed to resolve whether the
#             gap between observed and permuted holds up statistically.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
#
# 2026-08-31  ADDED gene-level permutation p-values (GTEx/FastQTL eGene
#             procedure) -> results/tables/eqtl_gene_level_perm.csv.
#
# 2026-08-31  CHANGED N_PERM default from 100 to 1000 (T21_N_PERM env var
#             still overrides). The planned Alpine offload that the 100
#             default was set for has been cancelled, and 1000 permutations
#             over the current 4-gene DE set is cheap locally (seconds),
#             restoring the previously committed 0.000999 p-value floor.
#
# 2026-08-31  RENAMED eqtl_gene_level_perm.csv column explained_perm ->
#             cis_eqtl_detected. Reason: the tight plan explicitly retires
#             any "explained by eQTL" claim (see docs/REPO_STATE.md decision
#             log); "explained" implies a causal/complete account this test
#             does not make. The column still means the same thing (gene-
#             level permutation q_gene_bh < FDR_GENE) under a name that
#             doesn't overclaim.
#
# 2026-09-01  CHANGED make_panel(): an empty gene set now removes any PDF
#             already at the declared output path instead of returning
#             silently, so an earlier run's panel cannot survive as a stale
#             output (CodeRabbit review, PR #3).
# 2026-09-10  REPLACED the negative controls of the retired any-variant rule
#             (direction flip, genotype shuffle; eqtl_negative_controls.csv)
#             with standalone controls of the classification test itself:
#             a negative control that pairs each deviating gene with the cis
#             variant set of a distant tested gene (decoy_min_distance) and a
#             positive control of strong GTEx eGenes from script 02, both run
#             through the same gene_level_tests() runner as the observed set
#             (scripts/lib/eqtl_controls.R). The observed test is unchanged:
#             same variants, seeds and BH, now called through the runner.
#             Positive-control rows are dropped from t21_dosage_per_variant.csv
#             and the representative-variant table so script 04 is unaffected.
# 2026-09-10  REMOVED the boxplot figures (t21_dosage_boxplots_*.pdf); the
#             per-gene dosage panels are drawn by scripts/11_eqtl_figures.R
#             from the best variant of the permutation test, with the control
#             sets alongside. Tables are unchanged.
