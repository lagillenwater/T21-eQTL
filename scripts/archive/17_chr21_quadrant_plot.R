# 17_chr21_quadrant_plot.R
#
# Purpose: Quadrant scatter for chr21 genes - ploidy-corrected deviation
#          (norm_log2FC) on x, within-T21 slope at the gene's
#          representative cis variant on y.
#
#          Quadrant interpretation:
#            Q1 (top-right):    pos deviation, pos T21 slope
#                               -> direction concordant (eQTL pushes the
#                                  same way the deviation goes)
#            Q3 (bottom-left):  neg deviation, neg T21 slope
#                               -> direction concordant
#            Q2 (top-left):     pos deviation, neg T21 slope
#                               -> direction discordant (top eQTL points
#                                  the opposite way)
#            Q4 (bottom-right): neg deviation, pos T21 slope
#                               -> direction discordant
#
#          By construction, supported genes always sit in Q1/Q3 (we pick
#          the representative by direction match); the value of this plot
#          is showing where the discordant genes (Q2/Q4) fall and how
#          tightly the supported genes cluster on the diagonal.
#
#          Within-T21 only (per user) - GTEx-slope panel removed.
#
# Inputs:
#   - results/tables/chr21_lane_assignments.csv
#   - results/tables/t21_dosage_per_variant.csv
#
# Outputs:
#   - results/figures/chr21_quadrant_plot.pdf
#   - results/figures/chr21_quadrant_plot.png
#   - results/figures/chr21_quadrant_plot_session_info.txt
#
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

set.seed(42)

cat("=== T21-eQTL: chr21 quadrant plot ===\n\n")

# =============================================================================
# STEP 1: Load lane assignments and per-variant fits
# =============================================================================

cat("Step 1: Loading lane assignments and per-variant fits...\n")

lanes <- fread("results/tables/chr21_lane_assignments.csv")
fits  <- fread("results/tables/t21_dosage_per_variant.csv")

# Pick representative variant per gene with this fallback chain:
#   1. strongest_supp_variant   (T21-reproducible direction match)
#   2. strongest_dir_variant    (GTEx direction match, no T21 reproduction)
#   3. strongest_overall_variant (smallest GTEx pval, regardless of dir)
lanes[, rep_variant := fcase(
  !is.na(strongest_supp_variant)    & nzchar(strongest_supp_variant),
    strongest_supp_variant,
  !is.na(strongest_dir_variant)     & nzchar(strongest_dir_variant),
    strongest_dir_variant,
  !is.na(strongest_overall_variant) & nzchar(strongest_overall_variant),
    strongest_overall_variant,
  default = NA_character_)]

# Pull the GTEx + within-T21 slopes for that representative
rep_fits <- merge(
  lanes[!is.na(rep_variant),
        .(ensembl_stable, rep_variant, sig_lane, eqtl_lane, norm_log2FC,
          deviation_vs_cohort_sd, Gene_name)],
  fits[, .(ensembl_stable, variant_id, gtex_slope, t21_slope, t21_p)],
  by.x = c("ensembl_stable", "rep_variant"),
  by.y = c("ensembl_stable", "variant_id"),
  all.x = TRUE)

cat(sprintf("  chr21 genes with a representative variant: %d\n",
            nrow(rep_fits)))

# Lane label for plotting (matches script 13)
rep_fits[, lane_label := fcase(
  eqtl_lane == "explained",   "eQTL-supported",
  eqtl_lane == "unexplained", "eQTL-tested, not supported",
  eqtl_lane == "no_GTEx_data", "no GTEx eQTL data")]
rep_fits[, lane_label := factor(lane_label,
  levels = c("eQTL-supported",
             "eQTL-tested, not supported",
             "no GTEx eQTL data"))]
rep_fits[, sig_class := fcase(
  sig_lane == "DE_low",  "DE_low",
  sig_lane == "DE_high", "DE_high",
  default               = "Not DE")]

# =============================================================================
# STEP 2: Build the two quadrant panels
# =============================================================================

cat("\nStep 2: Building quadrant panels...\n")

palette_lane <- c(
  "eQTL-supported"             = "#4DAF4A",
  "eQTL-tested, not supported" = "#FF7F0E",
  "no GTEx eQTL data"          = "#999999",
  "below cohort noise"         = "#CFD8DC")

# Bring in baseMean for the within-T21 rescale and add a magnitude flag
rep_fits <- merge(rep_fits,
                  lanes[, .(ensembl_stable, baseMean,
                            passes_magnitude_filter)],
                  by = "ensembl_stable")
rep_fits[, t21_slope_pct := t21_slope / (baseMean + 1)]

# Add the new "below cohort noise" lane label for genes that didn't make
# the magnitude filter
rep_fits[passes_magnitude_filter == FALSE,
         lane_label := "below cohort noise"]
rep_fits[, lane_label := factor(lane_label,
  levels = c("eQTL-supported",
             "eQTL-tested, not supported",
             "no GTEx eQTL data",
             "below cohort noise"))]

# Label DE genes that pass the magnitude filter (these are the ones the
# user actually cares about - the rest are intentionally greyed)
de_to_label <- rep_fits[sig_class != "Not DE" &
                        passes_magnitude_filter == TRUE &
                        !is.na(t21_slope_pct)]
cat(sprintf("  Genes flagged for labeling (DE, passes mag filter): %d\n",
            nrow(de_to_label)))

quadrant_layer <- function(p) {
  p +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey20") +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey20") +
    annotate("text", x =  0.05, y =  Inf, hjust = 0, vjust = 1.4,
             label = "Q1: concordant (+/+)", size = 2.7,
             color = "grey30") +
    annotate("text", x =  0.05, y = -Inf, hjust = 0, vjust = -0.6,
             label = "Q4: discordant (-T21 slope, +dev)", size = 2.7,
             color = "grey30") +
    annotate("text", x = -0.05, y =  Inf, hjust = 1, vjust = 1.4,
             label = "Q2: discordant (+T21 slope, -dev)", size = 2.7,
             color = "grey30") +
    annotate("text", x = -0.05, y = -Inf, hjust = 1, vjust = -0.6,
             label = "Q3: concordant (-/-)", size = 2.7,
             color = "grey30")
}

p_t21 <- ggplot(rep_fits[!is.na(t21_slope_pct)],
                aes(x = norm_log2FC, y = t21_slope_pct,
                    color = lane_label, shape = sig_class)) +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_text_repel(data = de_to_label,
                  aes(label = Gene_name), size = 2.5,
                  max.overlaps = 30, color = "grey20") +
  scale_color_manual(values = palette_lane) +
  scale_shape_manual(values = c("DE_low" = 16, "DE_high" = 17,
                                "Not DE" = 1)) +
  coord_cartesian(ylim = quantile(rep_fits$t21_slope_pct,
                                  c(0.02, 0.98), na.rm = TRUE)) +
  labs(
    title    = "Within-T21 slope vs ploidy-corrected deviation",
    subtitle = paste0("Per-gene representative cis variant. Y is the ",
                      "within-T21 raw-count slope rescaled by ",
                      "baseMean+1 for cross-gene comparison (2-98% ",
                      "range shown). Concordant quadrants (Q1, Q3) are ",
                      "where the within-T21 dosage effect points the ",
                      "same way as the deviation."),
    x        = "norm log2(T21 / Control)  [ploidy-corrected]",
    y        = "within-T21 slope / (baseMean + 1)  [per allele copy]",
    color    = NULL, shape = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        legend.box = "vertical")

p_t21 <- quadrant_layer(p_t21)

# =============================================================================
# STEP 3: Compose and save
# =============================================================================

cat("\nStep 3: Saving panel...\n")

ggsave("results/figures/chr21_quadrant_plot.pdf", p_t21,
       width = 11, height = 7, bg = "white")
ggsave("results/figures/chr21_quadrant_plot.png", p_t21,
       width = 11, height = 7, dpi = 150, bg = "white")
cat("  Saved chr21_quadrant_plot.{pdf,png}\n")

# =============================================================================
# STEP 4: Verification (per-quadrant gene counts)
# =============================================================================

cat("\n=== Verification ===\n")
rep_fits[, q_t21 := fcase(
  is.na(t21_slope_pct) | is.na(norm_log2FC), "no_data",
  norm_log2FC > 0 & t21_slope_pct > 0,       "Q1_concordant",
  norm_log2FC < 0 & t21_slope_pct < 0,       "Q3_concordant",
  norm_log2FC > 0 & t21_slope_pct < 0,       "Q2_discordant",
  norm_log2FC < 0 & t21_slope_pct > 0,       "Q4_discordant",
  default                                  = "axis")]

cat("\nWithin-T21 slope quadrant counts by lane (DE genes only):\n")
print(rep_fits[sig_class != "Not DE",
               .N, by = .(sig_class, lane_label, q_t21)
              ][order(sig_class, lane_label, q_t21)])

writeLines(capture.output(sessionInfo()),
           "results/figures/chr21_quadrant_plot_session_info.txt")

cat("\n=== Quadrant plot complete ===\n")
