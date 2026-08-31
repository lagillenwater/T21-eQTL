# 11_t21_dosage_boxplots.R
#
# Purpose: Within-T21 directional concordance + visualization. For each target
#          gene, regress expression on alt-allele dosage in T21 subjects only,
#          flag variants where the within-T21 slope sign matches the GTEx
#          slope sign ("supportive eQTL"), pick the representative variant
#          per gene as the most significant supportive eQTL, and faceted-
#          boxplot expression by genotype for genes that have at least one
#          supportive variant.
#
# Inputs:
#   - data/processed/genotypes_filtered.csv
#   - data/processed/eqtl_target_variants.csv
#   - data/processed/eqtl_supported_genes.csv
#   - data/processed/count_matrix.csv
#   - data/processed/sample_metadata.csv
#
# Outputs:
#   - results/tables/t21_dosage_per_variant.csv
#   - results/tables/t21_representative_variants.csv
#   - results/figures/t21_dosage_boxplots_de_low_fc.pdf
#   - results/figures/t21_dosage_boxplots_sig_high_fc.pdf
#   - results/figures/t21_dosage_boxplots_session_info.txt
#
# Date: 2026-04-30

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source("scripts/lib/eqtl_fit.R")

set.seed(42)

ALPHA_REPRO <- 0.05   # within-T21 nominal p for a cis variant to count as
                      # reproducible; must match scripts/04_chr21_lane_assignment.R

cat("=== T21-eQTL: Within-T21 Dosage-Expression Boxplots ===\n\n")

# =============================================================================
# STEP 1: Load inputs and restrict to T21
# =============================================================================

cat("Step 1: Loading inputs and restricting to T21...\n")

geno    <- fread("data/processed/genotypes_filtered.csv")
targets <- fread("data/processed/eqtl_target_variants.csv")
genes   <- fread("data/processed/eqtl_supported_genes.csv")
meta    <- fread("data/processed/sample_metadata.csv")
counts  <- fread("data/processed/count_matrix.csv")

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

counts[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
target_stable <- unique(sub("\\..*$", "", targets$gene_id))
counts_target <- counts[ensembl_stable %in% target_stable]
cat(sprintf("  Target genes matched in count matrix: %d / %d\n",
            uniqueN(counts_target$ensembl_stable), length(target_stable)))

t21_labids <- intersect(colnames(counts_target),
                        subj_to_labid[shared_subjects])

expr_long <- melt(counts_target,
                  id.vars = c("EnsemblID", "Gene_name", "Chr",
                              "ensembl_stable"),
                  measure.vars = t21_labids,
                  variable.name = "LabID",
                  value.name = "count",
                  variable.factor = FALSE)

labid_to_subj <- setNames(meta_t21$subject_id, meta_t21$LabID)
expr_long[, subject_id := labid_to_subj[LabID]]

cat(sprintf("  Genes with expression rows:   %d\n",
            uniqueN(expr_long$EnsemblID)))
cat(sprintf("  Subjects in expression long:  %d\n",
            uniqueN(expr_long$subject_id)))

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
  expr_long[, .(ensembl_stable, subject_id, count)],
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
  y <- count
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

cat(sprintf("  Variants tested: %d  Supportive: %d\n",
            nrow(fit_table), sum(fit_table$supportive, na.rm = TRUE)))

fwrite(fit_table, "results/tables/t21_dosage_per_variant.csv")
stopifnot(file.exists("results/tables/t21_dosage_per_variant.csv"))

# =============================================================================
# NEGATIVE CONTROLS for the "explained" call
# =============================================================================
# C2 direction flip       - negate the deviation direction and recount. The rule
#                           only asks that some cis variant point the same way
#                           as the deviation, so flipping measures its
#                           discriminating power.
# C1 genotype permutation - shuffle subject labels and refit. Destroys the
#                           genotype-expression link while preserving the number
#                           of variants per gene and their LD, which is the
#                           multiplicity the "any variant" rule ignores.
# A sound criterion collapses toward 0% under both.

cat("\n=== Negative controls ===\n")
N_PERM_SETS <- as.integer(Sys.getenv("T21_PERM_SETS", "20"))

score_explained <- function(dt, deviation_sign, t21_slope, t21_p) {
  supportive <- sign(dt$gtex_slope) == deviation_sign &
                sign(t21_slope) == sign(dt$gtex_slope) &
                !is.na(t21_p) & t21_p < ALPHA_REPRO
  tapply(supportive, dt$Gene_name, function(x) any(x, na.rm = TRUE))
}

obs  <- score_explained(fit_table, fit_table$observed_direction,
                        fit_table$t21_slope, fit_table$t21_p)
flip <- score_explained(fit_table, -fit_table$observed_direction,
                        fit_table$t21_slope, fit_table$t21_p)
cat(sprintf("  observed:          %d/%d explained (%.1f%%)\n",
            sum(obs), length(obs), 100 * mean(obs)))
cat(sprintf("  direction flipped: %d/%d explained (%.1f%%)\n",
            sum(flip), length(flip), 100 * mean(flip)))

gwide <- dcast(geno_t21, subject_id ~ variant_id, value.var = "alt_dosage")
subj  <- gwide$subject_id
G_all <- as.matrix(gwide[, -1])
de_genes <- unique(fit_table$Gene_name)

perm_rate <- vapply(seq_len(N_PERM_SETS), function(b) {
  set.seed(1000 + b)
  ord <- sample.int(length(subj))
  mean(vapply(de_genes, function(g) {
    vg   <- fit_table[Gene_name == g]
    cols <- intersect(vg$variant_id, colnames(G_all))
    if (length(cols) == 0) return(FALSE)
    e <- expr_of_gene(g, subj, counts, meta_t21)
    if (all(is.na(e))) return(FALSE)
    fits <- fit_variants(G_all[ord, cols, drop = FALSE], e)
    vgm  <- vg[match(cols, vg$variant_id)]
    any(sign(vgm$gtex_slope) == vgm$observed_direction &
        sign(fits$slope) == sign(vgm$gtex_slope) &
        !is.na(fits$p) & fits$p < ALPHA_REPRO, na.rm = TRUE)
  }, logical(1)))
}, numeric(1))
cat(sprintf("  genotype permuted: %.1f%% explained (mean of %d shuffles, range %.1f-%.1f%%)\n",
            100 * mean(perm_rate), N_PERM_SETS,
            100 * min(perm_rate), 100 * max(perm_rate)))

# n_explained for the genotype_permutation row is a MEAN across shuffles
# (mean(perm_rate) * n genes), not an integer gene count like the other two
# rows -- rounding it to an integer here previously produced a row like
# "0 of 3 = 15.02%", self-contradictory against its own pct_explained.
neg_ctrl <- data.table(
  control        = c("observed", "direction_flip", "genotype_permutation"),
  n_genes_tested = c(length(obs), length(flip), length(de_genes)),
  n_explained    = c(sum(obs), sum(flip), mean(perm_rate) * length(de_genes)),
  pct_explained  = c(100 * mean(obs), 100 * mean(flip), 100 * mean(perm_rate)))
fwrite(neg_ctrl, "results/tables/eqtl_negative_controls.csv")
stopifnot(file.exists("results/tables/eqtl_negative_controls.csv"))
cat("  Wrote results/tables/eqtl_negative_controls.csv\n")

# =============================================================================
# STEP 4: Pick representative variant per gene (most sig supportive eQTL)
# =============================================================================

cat("\nStep 4: Selecting representative variant per gene...\n")

supportive <- fit_table[supportive == TRUE & !is.na(t21_p)]
setorder(supportive, ensembl_stable, t21_p)
representatives <- supportive[, head(.SD, 1L), by = ensembl_stable]

fwrite(representatives,
       "results/tables/t21_representative_variants.csv")
stopifnot(file.exists("results/tables/t21_representative_variants.csv"))

n_rep_low  <- sum(representatives$gene_set == "DE_low_FC")
n_rep_high <- sum(representatives$gene_set == "Sig_high_FC")
cat(sprintf("  Representative-variant genes:  %d total ",
            nrow(representatives)))
cat(sprintf("(DE_low_FC: %d / 66, Sig_high_FC: %d / 16)\n",
            n_rep_low, n_rep_high))

# =============================================================================
# STEP 5: Boxplot data and figure
# =============================================================================

cat("\nStep 5: Building boxplot panels...\n")

plot_df <- merge(
  geno_expr[, .(variant_id, ensembl_stable, Gene_name, gene_set,
                subject_id, alt_dosage, count)],
  representatives[, .(variant_id, ensembl_stable)],
  by = c("variant_id", "ensembl_stable"), all = FALSE
)

plot_df[, dosage_factor := factor(alt_dosage, levels = 0:3)]

# Annotation: facet label = "GENE\nvariant_id\nGTEx slope=..., T21 p=..."
panel_labels <- representatives[, .(
  ensembl_stable, Gene_name, gene_set,
  variant_id, gtex_slope, t21_slope, t21_p, t21_n
)]
panel_labels[, facet_lab := sprintf(
  "%s\n%s\nGTEx %.2f | T21 b=%.1f p=%.1g",
  Gene_name, variant_id, gtex_slope, t21_slope, t21_p)]

plot_df <- merge(plot_df, panel_labels[, .(ensembl_stable, facet_lab)],
                 by = "ensembl_stable")

make_panel <- function(df, title, out_pdf) {
  if (nrow(df) == 0) {
    cat(sprintf("  No genes for '%s' - skipping.\n", title))
    return(invisible())
  }
  ng <- uniqueN(df$Gene_name)
  ncol <- min(5, ceiling(sqrt(ng)))
  nrow_pl <- ceiling(ng / ncol)

  p <- ggplot(df, aes(x = dosage_factor, y = count + 1)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey40") +
    geom_jitter(width = 0.18, size = 0.6, alpha = 0.4) +
    facet_wrap(~ facet_lab, scales = "free_y", ncol = ncol) +
    scale_y_continuous(trans = "log10") +
    labs(title = title,
         subtitle = "T21 subjects only; representative supportive eQTL per gene",
         x = "Alt-allele dosage in T21 (0-3)",
         y = "Raw expression count + 1 (log10)") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7, lineheight = 0.9),
          plot.title = element_text(face = "bold"))

  w <- max(8, ncol * 2.2)
  h <- max(5, nrow_pl * 2.2)
  ggsave(out_pdf, p, width = w, height = h, limitsize = FALSE)
  cat(sprintf("  Saved: %s (%d genes, %.1f x %.1f in)\n",
              out_pdf, ng, w, h))
}

make_panel(
  plot_df[gene_set == "DE_low_FC"],
  "DE Low FC genes - within-T21 dosage vs expression",
  "results/figures/t21_dosage_boxplots_de_low_fc.pdf"
)

make_panel(
  plot_df[gene_set == "Sig_high_FC"],
  "Sig >=1.5 FC genes - within-T21 dosage vs expression",
  "results/figures/t21_dosage_boxplots_sig_high_fc.pdf"
)

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
           "results/figures/t21_dosage_boxplots_session_info.txt")

cat("\n=== Boxplot script complete ===\n")

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
