# 20_combined_figures.R
#
# Purpose: The three main manuscript figures, each assembled from the panels
#          the source scripts save with save_panels() (scripts/lib/figure_style.R).
#          Nothing is recomputed here; rerun a source script to change a panel.
#
#   fig1_dosage_correction            A-B  volcano, uncorrected / ploidy-corrected (script 07)
#                                     C-D  log2FC density / ECDF, HSA21 and HSA22 (script 06)
#   fig2_classification_and_genotype  A    lane flow (script 12)
#                                     B    HSA21 map of the deviating genes (script 12)
#                                     C-D  expression by genotype, higher / lower than expected (script 11)
#   fig3_eqtl_effects_and_frequencies A-B  T21/GTEx effect ratio; spike-in power (script 18)
#                                     C-E  allele frequencies; deviation and reach; within-T21
#                                          variance explained (script 13); A-C on the top row
#
# Panel letters run continuously through each figure; the source panels' own
# tags are replaced. Sizes are for a full journal page width (7.5 in).
#
# Inputs (per run): figures/panels/{volcano_chr21, ploidy_correction_distributions,
#   chr21_flow_and_map, eqtl_dosage_panels, eqtl_detection_overview,
#   af_deviation_bound}.rds
# Outputs (per run): figures/fig1_dosage_correction.{pdf,png},
#   figures/fig2_classification_and_genotype.{pdf,png},
#   figures/fig3_eqtl_effects_and_frequencies.{pdf,png}
#
# Usage: Rscript scripts/20_combined_figures.R --run adjusted
#   after scripts 06, 07, 11, 12, 13 and 18 for the same run.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggsankey)
  library(patchwork)
})
source("scripts/lib/run.R"); run <- load_run()

cat(sprintf("=== T21-eQTL: combined main figures [run: %s] ===\n\n", run$name))

panels <- function(stem) {
  path <- file.path(run$dirs$figures, "panels", paste0(stem, ".rds"))
  if (!file.exists(path)) stop("missing ", path, ": run the script that draws ", stem, " first", call. = FALSE)
  readRDS(path)
}

# One panel-letter style for every figure. Tags sit in the plot margin, so a
# collected axis title (Figure 2C-D) does not run under them.
TAG_THEME <- theme(plot.tag = element_text(face = "bold", size = 11), plot.tag.location = "margin")

save_fig <- function(fig, stem, width, height) {
  out <- run$figure(stem)
  ggsave(paste0(out, ".pdf"), fig, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(out, ".png"), fig, width = width, height = height, dpi = 300, limitsize = FALSE)
  for (f in paste0(out, c(".pdf", ".png"))) {
    if (!file.exists(f)) stop("failed to write ", f)
    cat(sprintf("  Saved: %s\n", f))
  }
}

# =============================================================================
# Figure 1: the ploidy correction
# =============================================================================
volc <- panels("volcano_chr21")
dist <- panels("ploidy_correction_distributions")

# One legend per row. The volcano legend sits in a guide area under the two
# volcanoes (as script 07 draws it): setting legend.position through `&`
# would switch panel B's legend back on and draw it twice.
# Both keys on one line.
volc_keys <- guides(
  colour = guide_legend(nrow = 1, order = 1, override.aes = list(size = 2.2, alpha = 1)),
  shape  = guide_legend(nrow = 1, order = 2, override.aes = list(size = 2.2, alpha = 1, colour = "grey30")))
row_volc <- (volc$uncorrected + volc_keys + theme(legend.position = "bottom")) +
  (volc$corrected + theme(legend.position = "none")) + guide_area() +
  plot_layout(design = "AB\nCC", guides = "collect", heights = c(1, 0.08)) &
  theme(legend.box = "horizontal", legend.box.just = "center", legend.spacing.x = unit(0.8, "lines"),
        legend.text = element_text(size = 7))
# Script 06 draws at base size 11; match the volcanoes' 9.
row_dist <- (dist$density + dist$ecdf + plot_layout(guides = "collect")) &
  theme(text = element_text(size = 9), legend.position = "bottom")

fig1 <- (row_volc / row_dist) + plot_layout(heights = c(1.2, 0.85)) +
  plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(fig1, "fig1_dosage_correction", width = 7.5, height = 7.5)

# =============================================================================
# Figure 2: classification of the HSA21 genes and expression by genotype
# =============================================================================
fm  <- panels("chr21_flow_and_map")
dos <- panels("eqtl_dosage_panels")

# Each block's own y title runs longer than the block at this size, and a
# shared one runs under the D tag: a short two-line title on each block (the
# caption says the covariate effects are removed).
Y_DOS <- expression(atop("Expression", "(" * log[2] * "-CPM, adjusted)"))
row_dos <- (dos$higher + labs(y = Y_DOS)) / (dos$lower + labs(y = Y_DOS))
# B keeps only the eQTL-outcome key; the caption explains the red and blue
# label colours (direction of deviation).
fig2 <- (fm$flow / (fm$map + guides(colour = "none")) / row_dos) +
  plot_layout(heights = c(1.5, 1.0, 1.85)) +
  plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(fig2, "fig2_classification_and_genotype", width = 7.5, height = 9)

# =============================================================================
# Figure 3: cis-eQTL effects in T21 and what they could do to the group ratio
# =============================================================================
eff <- panels("eqtl_detection_overview")
afb <- panels("af_deviation_bound")

# Top row: A, B and C (the allele frequencies). A's legend, on one line,
# serves the row (direction colours with the Expected-dosage grey, and the
# detection shapes); B repeats them and C has none.
eff_keys <- guides(colour = guide_legend(nrow = 1, order = 1),
                   shape  = guide_legend(nrow = 1, order = 2))
# Each row's keys sit in a guide area across the whole row, so a one-line
# legend is centred on the page rather than on the panels.
LEGEND_LINE <- theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 7),
                     legend.key.spacing.x = unit(2, "pt"), legend.spacing.x = unit(6, "pt"))
# Tag styles are set per panel here, not with `& TAG_THEME`, so E can differ.
# E has no y-axis text, so a margin tag would land left of its panel, over D's
# plot area: E's tag sits over its own panel's top-left corner instead.
TAG_E <- theme(plot.tag = element_text(face = "bold", size = 11, hjust = 0.5, vjust = -0.5),
               plot.tag.location = "panel", plot.tag.position = c(0, 1))
row_top <- ((eff$ratio + eff_keys + TAG_THEME) +
              (eff$power + guides(colour = "none", shape = "none") + TAG_THEME) + (afb$af + TAG_THEME) +
              guide_area() +
              plot_layout(design = "ABC\nDDD", widths = c(1, 1.1, 0.95), heights = c(1, 0.07),
                          guides = "collect")) & LEGEND_LINE
# Bottom row: D and E. D's keys on one line; its direction colours are the
# top row's.
bound_keys <- guides(colour = "none",
                     shape = guide_legend(nrow = 1, order = 1,
                                          override.aes = list(colour = c("grey20", "black"), fill = c(NA, "black"))),
                     linewidth = guide_legend(order = 2))
row_bottom <- ((afb$bound + bound_keys + TAG_THEME) + (afb$r2 + TAG_E) + guide_area() +
                 plot_layout(design = "AB\nCC", widths = c(2.4, 1), heights = c(1, 0.06), guides = "collect")) &
  LEGEND_LINE

# free(): D's long row labels would otherwise set the left edge of A's panel.
fig3 <- (row_top / free(row_bottom)) + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "A")
save_fig(fig3, "fig3_eqtl_effects_and_frequencies", width = 7.5, height = 7.75)

writeLines(capture.output(sessionInfo()), run$figure("combined_figures_session_info.txt"))
cat("\n=== Combined figures complete ===\n")
