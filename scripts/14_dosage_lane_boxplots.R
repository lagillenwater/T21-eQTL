# 14_dosage_lane_boxplots.R
#
# Purpose: Focused within-T21 dosage-vs-expression boxplots, one PDF per
#          non-empty (sig_lane, eqtl_lane) quadrant from
#          chr21_lane_assignments.csv. For explained genes, the panel uses
#          the strongest supportive-with-reproducibility variant per gene.
#          For unexplained genes, the panel uses the strongest GTEx direction-
#          matching variant (so the figure shows the "best try" that failed
#          to reproduce in T21).
#
#          Each panel is annotated with:
#            n_repro / n_dir_match / n_cis_total
#            norm_log2FC   (observed deviation from 1.5x ploidy)
#            GTEx slope and pval at the displayed variant
#            within-T21 slope and pval at the displayed variant
#
# Inputs:
#   - results/tables/chr21_lane_assignments.csv
#   - results/tables/t21_dosage_per_variant.csv
#   - data/processed/genotypes_filtered.csv
#   - data/processed/count_matrix.csv
#   - data/processed/sample_metadata.csv
#
# Outputs:
#   - results/figures/chr21_dosage_de_low_supported.pdf
#   - results/figures/chr21_dosage_de_low_not_supported.pdf
#   - results/figures/chr21_dosage_de_high_supported.pdf
#   - results/figures/chr21_dosage_de_high_not_supported.pdf  (if non-empty)
#   - results/figures/chr21_dosage_lane_boxplots_session_info.txt
#
# Author: Claude Code
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)

cat("=== T21-eQTL: Dosage Boxplots by Lane ===\n\n")

# =============================================================================
# STEP 1: Load lane assignments and pick representative variants per gene
# =============================================================================

cat("Step 1: Loading lane assignments...\n")

lanes  <- fread("results/tables/chr21_lane_assignments.csv")
fits   <- fread("results/tables/t21_dosage_per_variant.csv")

# For each gene-by-quadrant, choose the representative variant:
#   explained  -> strongest_supp_variant   (within-T21 reproducible)
#   unexplained -> strongest_dir_variant   (GTEx dir match; failed T21 repro)
# Genes in the no_GTEx_data lane have nothing to plot - skipped here.
panel_specs <- lanes[
  eqtl_lane %in% c("explained", "unexplained") &
    sig_lane %in% c("DE_low", "DE_high"),
  .(ensembl_stable, Gene_name, baseMean, raw_log2FC, norm_log2FC, norm_padj,
    sig_lane, eqtl_lane,
    n_cis_total, n_dir_match, n_supp_with_repro,
    representative_variant = fcase(
      eqtl_lane == "explained",                  strongest_supp_variant,
      eqtl_lane == "unexplained" &
        !is.na(strongest_dir_variant) &
        nzchar(strongest_dir_variant),           strongest_dir_variant,
      eqtl_lane == "unexplained",                strongest_overall_variant))]

# Sanity check: we should have a representative for every plotable gene.
n_missing <- sum(is.na(panel_specs$representative_variant))
if (n_missing > 0L) {
  cat(sprintf("  WARNING: %d genes lack a representative variant; skipping\n",
              n_missing))
  panel_specs <- panel_specs[!is.na(representative_variant)]
}

cat(sprintf("  Genes to plot: %d\n", nrow(panel_specs)))
print(panel_specs[, .N, by = .(sig_lane, eqtl_lane)])

# =============================================================================
# STEP 2: Build long-format expression x dosage at representative variants
# =============================================================================

cat("\nStep 2: Loading expression and genotypes...\n")

meta   <- fread("data/processed/sample_metadata.csv")
counts <- fread("data/processed/count_matrix.csv")
geno   <- fread("data/processed/genotypes_filtered.csv")

meta_t21 <- meta[Karyotype == "T21"]
labid_to_subj <- setNames(sub("[A-Z][0-9]*$", "", meta_t21$LabID),
                          meta_t21$LabID)
subj_to_labid <- setNames(meta_t21$LabID,
                          sub("[A-Z][0-9]*$", "", meta_t21$LabID))

# T21 genotypes only, restricted to representative variants
geno_t21 <- geno[karyotype == "T21" & !is.na(alt_dosage) &
                 variant_id %in% panel_specs$representative_variant]
geno_t21 <- geno_t21[, .(alt_dosage = mean(alt_dosage)),
                     by = .(variant_id, subject_id)]
cat(sprintf("  T21 subjects with genotype at any rep variant: %d\n",
            uniqueN(geno_t21$subject_id)))

# Expression in T21 for the genes we want to plot
counts[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
target_stable <- unique(panel_specs$ensembl_stable)
counts_target <- counts[ensembl_stable %in% target_stable]

t21_labids <- intersect(colnames(counts_target),
                        unname(subj_to_labid))
expr_long <- melt(counts_target,
                  id.vars = c("EnsemblID", "Gene_name", "Chr",
                              "ensembl_stable"),
                  measure.vars = t21_labids,
                  variable.name = "LabID",
                  value.name = "count",
                  variable.factor = FALSE)
expr_long[, subject_id := labid_to_subj[LabID]]

# Join genotype (per representative variant) to expression (per gene)
plot_long <- merge(
  panel_specs[, .(ensembl_stable, Gene_name_lane = Gene_name,
                  representative_variant, sig_lane, eqtl_lane,
                  baseMean, raw_log2FC, norm_log2FC, norm_padj,
                  n_cis_total, n_dir_match, n_supp_with_repro)],
  expr_long[, .(ensembl_stable, subject_id, count)],
  by = "ensembl_stable", allow.cartesian = TRUE
)
plot_long <- merge(
  plot_long,
  geno_t21,
  by.x = c("representative_variant", "subject_id"),
  by.y = c("variant_id",             "subject_id"),
  all = FALSE
)

# Per-variant within-T21 fit stats for the panel annotation
fit_keys <- unique(panel_specs[, .(ensembl_stable,
                                   variant_id = representative_variant)])
fit_anno <- merge(fit_keys,
                  fits[, .(ensembl_stable, variant_id, gtex_slope, gtex_pval,
                           t21_slope, t21_p, t21_n)],
                  by = c("ensembl_stable", "variant_id"))

panel_specs <- merge(panel_specs, fit_anno,
                     by.x = c("ensembl_stable", "representative_variant"),
                     by.y = c("ensembl_stable", "variant_id"),
                     all.x = TRUE)

# Build facet label per gene
panel_specs[, facet_lab := sprintf(
  paste0("%s | norm_log2FC=%.2f\n%s\n",
         "n_repro=%d / n_dir=%d / n_cis=%d\n",
         "GTEx b=%.2f p=%.1g | T21 b=%.1f p=%.1g"),
  Gene_name, norm_log2FC, representative_variant,
  n_supp_with_repro, n_dir_match, n_cis_total,
  gtex_slope, gtex_pval, t21_slope, t21_p)]

plot_long <- merge(plot_long,
                   panel_specs[, .(ensembl_stable, facet_lab)],
                   by = "ensembl_stable")
plot_long[, dosage_factor := factor(alt_dosage, levels = 0:3)]

# =============================================================================
# STEP 3: Render one PDF per (sig_lane, eqtl_lane) quadrant
# =============================================================================

cat("\nStep 3: Rendering per-quadrant PDFs...\n")

render_panel <- function(df, out_pdf, title, subtitle) {
  if (nrow(df) == 0L) {
    cat(sprintf("  Skipping (empty): %s\n", out_pdf))
    return(invisible())
  }
  ng   <- uniqueN(df$ensembl_stable)
  ncol <- min(5L, ceiling(sqrt(ng)))
  nrow_pl <- ceiling(ng / ncol)

  p <- ggplot(df, aes(x = dosage_factor, y = count + 1)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey40") +
    geom_jitter(width = 0.18, size = 0.5, alpha = 0.4) +
    facet_wrap(~ facet_lab, scales = "free_y", ncol = ncol) +
    scale_y_continuous(trans = "log10") +
    labs(title    = title,
         subtitle = subtitle,
         x = "Alt-allele dosage in T21 (0-3)",
         y = "Raw expression count + 1 (log10)") +
    theme_bw(base_size = 8) +
    theme(strip.text = element_text(size = 6, lineheight = 0.9),
          plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9))

  w <- max(8, ncol * 2.4)
  h <- max(5, nrow_pl * 2.4)
  ggsave(out_pdf, p, width = w, height = h, limitsize = FALSE)
  cat(sprintf("  Saved %s (%d genes, %.1f x %.1f in)\n",
              out_pdf, ng, w, h))
}

quadrants <- list(
  list(sig = "DE_low",  eqtl = "explained",
       file = "results/figures/chr21_dosage_de_low_supported.pdf",
       title = paste0("DE_low - eQTL-SUPPORTED (representative = strongest ",
                      "supportive-with-T21-reproducibility variant)"),
       sub   = paste0("Within-T21 alt-allele carriers have lower expression",
                      "; cohort sits below 1.5x ploidy ceiling")),
  list(sig = "DE_low",  eqtl = "unexplained",
       file = "results/figures/chr21_dosage_de_low_not_supported.pdf",
       title = paste0("DE_low - eQTL-TESTED, NOT SUPPORTED (representative ",
                      "= strongest GTEx direction match; failed T21 ",
                      "reproducibility)"),
       sub   = paste0("Tested cis variants do not reproduce in T21; ",
                      "interpret with deviation_vs_cohort_sd - small ",
                      "magnitude here is consistent with cohort noise, ",
                      "not necessarily evidence of compensation")),
  list(sig = "DE_high", eqtl = "explained",
       file = "results/figures/chr21_dosage_de_high_supported.pdf",
       title = paste0("DE_high - eQTL-SUPPORTED (representative = ",
                      "strongest supportive-with-T21-reproducibility ",
                      "variant)"),
       sub   = paste0("Within-T21 alt-allele carriers have higher ",
                      "expression; cohort sits above 1.5x ploidy ceiling")),
  list(sig = "DE_high", eqtl = "unexplained",
       file = "results/figures/chr21_dosage_de_high_not_supported.pdf",
       title = paste0("DE_high - eQTL-TESTED, NOT SUPPORTED"),
       sub   = paste0("Tested cis variants do not reproduce in T21; ",
                      "interpret with deviation_vs_cohort_sd"))
)

for (q in quadrants) {
  sub_df <- plot_long[
    ensembl_stable %in% panel_specs[sig_lane == q$sig &
                                    eqtl_lane == q$eqtl, ensembl_stable]]
  render_panel(sub_df, q$file, q$title, q$sub)
}

# =============================================================================
# STEP 4: Verification
# =============================================================================

cat("\n=== Verification ===\n")
verify <- panel_specs[, .(n_genes = .N), by = .(sig_lane, eqtl_lane)]
print(verify[order(sig_lane, eqtl_lane)])

writeLines(capture.output(sessionInfo()),
           "results/figures/chr21_dosage_lane_boxplots_session_info.txt")

cat("\n=== Lane boxplots complete ===\n")
