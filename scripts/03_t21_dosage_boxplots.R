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

set.seed(42)

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

# Vectorized per-(variant, gene) linear regression via formulas: with one
# predictor we can compute slope, intercept, t, and p in closed form.
fit_table <- geno_expr[, {
  x  <- alt_dosage
  y  <- count
  n  <- .N
  if (n < 5 || var(x) == 0) {
    .(t21_n = n, t21_slope = NA_real_, t21_se = NA_real_,
      t21_t = NA_real_, t21_p = NA_real_)
  } else {
    mx   <- mean(x); my <- mean(y)
    sxx  <- sum((x - mx)^2)
    sxy  <- sum((x - mx) * (y - my))
    b    <- sxy / sxx
    a    <- my - b * mx
    res  <- y - (a + b * x)
    sse  <- sum(res^2)
    sigma2 <- sse / (n - 2)
    se_b <- sqrt(sigma2 / sxx)
    tval <- b / se_b
    pval <- 2 * pt(-abs(tval), df = n - 2)
    .(t21_n = n, t21_slope = b, t21_se = se_b,
      t21_t = tval, t21_p = pval)
  }
}, by = .(variant_id, ensembl_stable, Gene_name, gene_set,
          gtex_slope, gtex_pval, observed_direction)]

fit_table[, supportive := !is.na(t21_slope) &
                          sign(t21_slope) == sign(gtex_slope) &
                          sign(t21_slope) != 0]

cat(sprintf("  Variants tested: %d  Supportive: %d\n",
            nrow(fit_table), sum(fit_table$supportive, na.rm = TRUE)))

fwrite(fit_table, "results/tables/t21_dosage_per_variant.csv")

# =============================================================================
# STEP 4: Pick representative variant per gene (most sig supportive eQTL)
# =============================================================================

cat("\nStep 4: Selecting representative variant per gene...\n")

supportive <- fit_table[supportive == TRUE & !is.na(t21_p)]
setorder(supportive, ensembl_stable, t21_p)
representatives <- supportive[, head(.SD, 1L), by = ensembl_stable]

fwrite(representatives,
       "results/tables/t21_representative_variants.csv")

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
