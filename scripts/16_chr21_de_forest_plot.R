# 16_chr21_de_forest_plot.R
#
# Purpose: Per-gene forest plot for chr21 DE genes (DE_low and DE_high).
#          Each gene shows its ploidy-normalized log2FC point estimate and
#          95% CI (estimate +/- 1.96 * lfcSE). Colored by eqtl_lane.
#          Lets the reader judge whether a deviation is well-estimated or
#          noisy: a tight CI excluding zero by a wide margin is a real
#          deviation; a CI that just barely excludes zero is borderline
#          regardless of padj.
#
#          Vertical reference at log2FC = 0 (the ploidy-corrected null) and
#          at +/- 1 cohort-noise SD lines (computed from non-chr21 genes)
#          to anchor magnitude vs typical cohort variation.
#
# Inputs:
#   - results/tables/chr21_lane_assignments.csv
#   - results/tables/deseq2_all_genes_ploidy_normalized.csv  (for lfcSE)
#
# Outputs:
#   - results/figures/chr21_de_forest_plot.pdf
#   - results/figures/chr21_de_forest_plot.png
#   - results/figures/chr21_de_forest_plot_session_info.txt
#
# Author: Claude Code
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)

cat("=== T21-eQTL: Chr21 DE forest plot ===\n\n")

# =============================================================================
# STEP 1: Load and join lane + lfcSE
# =============================================================================

cat("Step 1: Loading lane assignments and lfcSE...\n")

lanes <- fread("results/tables/chr21_lane_assignments.csv")
all_lfc <- fread("results/tables/deseq2_all_genes_ploidy_normalized.csv")
all_lfc[, ensembl_stable := sub("\\..*$", "", EnsemblID)]

m <- merge(lanes,
           all_lfc[, .(ensembl_stable, lfcSE)],
           by = "ensembl_stable", all.x = TRUE)

# Restrict to DE chr21 genes that pass paper's filters
de <- m[sig_lane %in% c("DE_low", "DE_high") &
        low_expr == FALSE & high_repeat == FALSE &
        !is.na(lfcSE)]
cat(sprintf("  DE chr21 genes after paper filters: %d\n", nrow(de)))

# Cohort-noise SD reference (matches what script 12 uses)
non_chr21_sd <- sd(
  all_lfc[Chr != "chr21" & !is.na(log2FoldChange), log2FoldChange])
cat(sprintf("  Cohort-noise SD (non-chr21): %.3f\n", non_chr21_sd))

# =============================================================================
# STEP 2: Build forest data
# =============================================================================

de[, lo95 := norm_log2FC - 1.96 * lfcSE]
de[, hi95 := norm_log2FC + 1.96 * lfcSE]

# Per user request: coloring keyed to whether the gene survives the
# within-cohort SD filter, not the eQTL lane. Genes that fail the filter
# are greyed out so the eye locks on the magnitude-meaningful set.
de[, color_group := fifelse(passes_magnitude_filter,
                            "passes magnitude filter",
                            "below cohort noise (greyed)")]
de[, color_group := factor(color_group,
  levels = c("passes magnitude filter",
             "below cohort noise (greyed)"))]
de[, sig_lane := factor(sig_lane, levels = c("DE_high", "DE_low"))]

# Order genes within each sig_lane by norm_log2FC
de[, gene_order := factor(Gene_name,
                          levels = de[order(sig_lane, norm_log2FC),
                                      Gene_name])]

cat(sprintf("\n  Genes plotted: %d (DE_low: %d, DE_high: %d)\n",
            nrow(de), sum(de$sig_lane == "DE_low"),
            sum(de$sig_lane == "DE_high")))

# =============================================================================
# STEP 3: Render
# =============================================================================

cat("\nStep 3: Rendering forest plot...\n")

p <- ggplot(de, aes(x = norm_log2FC, y = gene_order, color = color_group)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey20") +
  geom_vline(xintercept = c(-non_chr21_sd, non_chr21_sd),
             linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = c(-2 * non_chr21_sd, 2 * non_chr21_sd),
             linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                width = 0, linewidth = 0.4, alpha = 0.7,
                orientation = "y") +
  geom_point(size = 1.4) +
  scale_color_manual(values = c(
    "passes magnitude filter"     = "#1F77B4",
    "below cohort noise (greyed)" = "grey75")) +
  facet_grid(sig_lane ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "chr21 DE genes: ploidy-normalized log2FC with 95% CI",
    subtitle = sprintf(paste0("padj < 0.01; CIs are estimate +/- ",
                              "1.96*lfcSE. Dashed lines = +/-1 cohort-",
                              "noise SD (%.2f); dotted = +/-2 SD. ",
                              "Blue = passes magnitude filter (eligible ",
                              "for eQTL testing); grey = below cohort ",
                              "noise (not eQTL-tested)."),
                       non_chr21_sd),
    x        = "norm log2(T21 / Control), ploidy-corrected",
    y        = NULL,
    color    = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position    = "top",
        axis.text.y        = element_text(size = 6),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

# Sizing: vertical scales with gene count
n_total <- nrow(de)
height <- max(8, 0.13 * n_total)

ggsave("results/figures/chr21_de_forest_plot.pdf", p,
       width = 10, height = height, bg = "white", limitsize = FALSE)
ggsave("results/figures/chr21_de_forest_plot.png", p,
       width = 10, height = height, dpi = 150, bg = "white",
       limitsize = FALSE)
cat(sprintf("  Saved chr21_de_forest_plot.{pdf,png} (%.1f in tall)\n",
            height))

# =============================================================================
# STEP 4: Verification
# =============================================================================

cat("\n=== Verification ===\n")
cat(sprintf("\nMagnitude-filter breakdown:\n"))
print(de[, .(n_genes = .N,
             median_lfcSE = round(median(lfcSE), 3),
             median_devmag = round(median(abs(norm_log2FC)), 3)),
         by = .(sig_lane, color_group)])

writeLines(capture.output(sessionInfo()),
           "results/figures/chr21_de_forest_plot_session_info.txt")

cat("\n=== Forest plot complete ===\n")
