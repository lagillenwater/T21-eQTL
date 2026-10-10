# 06_chr21_distribution_panel.R
#
# Purpose: Show what the ploidy correction does at the population level, and
#          that it acts only where it should. Two panels (density, ECDF) of
#          the T21-vs-Control log2 fold change for chr21 genes of the target
#          biotypes (scripts/lib/biotypes.R)
#          on both scales - uncorrected, where the extra copy puts the
#          distribution at log2(1.5), and ploidy-corrected, where a gene that
#          follows dosage expectation sits at 0 - with chr22 on both scales as
#          the control. The correction only rescales chr21 counts, so the two
#          chr22 curves must coincide.
#
#          Scope what this shows: it is the population-level backdrop for the
#          per-gene classification (script 04). A corrected chr21 curve that
#          overlaps the chr22 curve is consistent with Hunter et al.'s "no
#          dosage compensation" at population level; per-gene compensation
#          candidates are not read from this figure.
#
# Inputs:
#   - results/tables/deseq2_all_genes_both_analyses.csv   (script 01)
#
# Outputs:
#   - results/figures/ploidy_correction_distributions.pdf
#   - results/figures/ploidy_correction_distributions.png
#   - results/tables/ploidy_correction_distribution_stats.csv
#   - results/figures/ploidy_correction_distributions_session_info.txt

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("scripts/lib/figure_style.R")   # strip_titles, record_titles
source("scripts/lib/ploidy_distributions.R")
source("scripts/lib/run.R"); run <- load_run()

# Match the rest of the pipeline (scripts 02 and 04): the chr21 gene set is
# TARGET_BIOTYPES, so the control chromosome is restricted the same way.
source("scripts/lib/biotypes.R")
CONTROL_CHR <- "chr22"
CHROMOSOMES <- c("chr21", CONTROL_CHR)
LFC_PLOIDY  <- log2(run$ploidy)
OUT_STEM <- run$figure("ploidy_correction_distributions")

# Chromosome is colour (the pair passes the colour-vision separation checks,
# CVD delta E >= 21 for every deficiency type). Scale is line weight: the
# uncorrected curve is a wide translucent band and the corrected curve a thin
# line drawn on top, so coincident curves (chr22) show as a line inside its
# band instead of one hiding the other.
# chr21 purple as in the all-genes volcano; the comparison chromosome grey
# (red and blue are reserved for the direction of deviation).
COL_CHR     <- setNames(c("#762A83", "grey40"), CHROMOSOMES)
SCALE_ORDER <- c("uncorrected", "ploidy-corrected")
LWD_SCALE   <- c("uncorrected" = 2.4, "ploidy-corrected" = 0.8)
ALPHA_SCALE <- c("uncorrected" = 0.35, "ploidy-corrected" = 1)

cat("=== T21-eQTL: Ploidy-correction distribution panel ===\n\n")

# =============================================================================
# STEP 1: Load both scales and build the long table
# =============================================================================

cat("Step 1: Loading DESeq2 results (both scales)...\n")
res  <- fread(run$table("deseq2_all_genes_both_analyses.csv"))
long <- ploidy_distribution_long(res, chromosomes = CHROMOSOMES,
                                 biotypes = TARGET_BIOTYPES)
n_chr <- long[scale == "uncorrected", .N, by = chromosome]
cat(sprintf("  %s: %d genes with estimates on both scales\n",
            n_chr$chromosome, n_chr$N), sep = "")

# =============================================================================
# STEP 2: Summary statistics
# =============================================================================

cat("\nStep 2: Summary statistics...\n")
stats <- ploidy_distribution_stats(long)
shift <- ploidy_shift_stats(long)
cat("\nPer chromosome x scale:\n"); print(stats)
cat("\nPer-gene shift (corrected - uncorrected):\n"); print(shift)

ctrl_shift <- max(abs(shift[chromosome == CONTROL_CHR, c(min_shift, max_shift)]))
if (ctrl_shift > 1e-3) {
  warning(sprintf(paste0("%s log2FC moves by up to %.3g under the ploidy ",
                         "correction; the control is expected to be unchanged"),
                  CONTROL_CHR, ctrl_shift))
}

lfc_of <- function(chr, sc) long[chromosome == chr & scale == sc, log2FC]
ks_corr <- ks.test(lfc_of("chr21", "ploidy-corrected"),
                   lfc_of(CONTROL_CHR, "ploidy-corrected"))
ks_raw  <- ks.test(lfc_of("chr21", "uncorrected"),
                   lfc_of(CONTROL_CHR, "uncorrected"))
cat(sprintf("\nKS chr21 vs %s, uncorrected:      D = %.3f  p = %.3g\n",
            CONTROL_CHR, ks_raw$statistic, ks_raw$p.value))
cat(sprintf("KS chr21 vs %s, ploidy-corrected: D = %.3f  p = %.3g\n",
            CONTROL_CHR, ks_corr$statistic, ks_corr$p.value))

stats_out <- merge(stats,
                   shift[, .(chromosome, median_shift, min_shift, max_shift)],
                   by = "chromosome")
setorder(stats_out, chromosome, scale)
fwrite(stats_out, run$table("ploidy_correction_distribution_stats.csv"))
stopifnot(file.exists(run$table("ploidy_correction_distribution_stats.csv")))
cat("  Wrote", run$table("ploidy_correction_distribution_stats.csv"), "\n")

# =============================================================================
# STEP 3: Density + ECDF
# =============================================================================

cat("\nStep 3: Building plots...\n")

x_lim <- c(min(-1.5, quantile(long$log2FC, 0.005)),
           max(1.5, quantile(long$log2FC, 0.995)))

style <- function(p, label_y = Inf, label_vjust = 1.6) {
  p +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = LFC_PLOIDY, colour = "grey60", linewidth = 0.4) +
    annotate("text", x = LFC_PLOIDY, y = label_y, label = "log[2]*(1.5)", parse = TRUE,
             hjust = -0.1, vjust = label_vjust, size = 3, colour = "grey30") +
    # Data keep the chr21/chr22 keys; the legend uses the manuscript's names.
    scale_colour_manual(values = COL_CHR, labels = c(chr21 = "HSA21", chr22 = "HSA22"),
                        name = "Chromosome") +
    scale_linewidth_manual(values = LWD_SCALE, breaks = SCALE_ORDER,
                           name = "Scale") +
    scale_alpha_manual(values = ALPHA_SCALE, breaks = SCALE_ORDER,
                       name = "Scale") +
    coord_cartesian(xlim = x_lim) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          legend.key.width = unit(1.6, "lines"))
}

# Uncorrected bands first, corrected lines on top (see the constants above).
corrected   <- long[scale == "ploidy-corrected"]
uncorrected <- long[scale == "uncorrected"]
aes_curve   <- aes(x = log2FC, colour = chromosome, linewidth = scale,
                   alpha = scale)

p_density <- style(
  ggplot(long, aes_curve) +
    # geom_line(stat = "density"), not geom_density: a density polygon applies
    # alpha to its fill, so the outline band would come out opaque.
    geom_line(stat = "density", data = uncorrected) +
    geom_line(stat = "density", data = corrected) +
    labs(title = "Density", x = expression(log[2]~"(T21 / D21)"), y = "Density")
)

p_ecdf <- style(
  ggplot(long, aes_curve) +
    stat_ecdf(data = uncorrected, geom = "step", pad = FALSE) +
    stat_ecdf(data = corrected,   geom = "step", pad = FALSE) +
    labs(title = "Empirical CDF", x = expression(log[2]~"(T21 / D21)"),
         y = "Cumulative fraction"),
  label_y = -Inf, label_vjust = -0.6
)

subtitle <- sprintf(
  paste0("Genes (", TARGET_BIOTYPES_LABEL, "): chr21 n = %d, %s n = %d. Median per-gene shift ",
         "from the correction: chr21 %.3f (log2(1.5) = %.3f), %s %.1e. ",
         "KS chr21 vs %s after correction: D = %.3f, p = %.2g"),
  n_chr[chromosome == "chr21", N], CONTROL_CHR, n_chr[chromosome == CONTROL_CHR, N],
  shift[chromosome == "chr21", median_shift], LFC_PLOIDY,
  CONTROL_CHR, shift[chromosome == CONTROL_CHR, median_shift],
  CONTROL_CHR, ks_corr$statistic, ks_corr$p.value)

# Panel letters only: the titles and subtitle are recorded for the captions.
TITLES <- file.path(run$dirs$figures, "removed_titles.txt")
record_titles("ploidy_correction_distributions", "",
              paste0("Ploidy correction shifts chr21 by log2(1.5) and leaves ", CONTROL_CHR, " unchanged"),
              subtitle, TITLES)
p_density <- strip_titles(p_density, "ploidy_correction_distributions", "A", TITLES, tag = TRUE)
p_ecdf    <- strip_titles(p_ecdf, "ploidy_correction_distributions", "B", TITLES, tag = TRUE)
# The two panels again, for the combined Figure 1 (script 20).
save_panels(list(density = p_density, ecdf = p_ecdf), "ploidy_correction_distributions", run)
panel <- (p_density | p_ecdf) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(OUT_STEM, ".pdf"), panel, width = 12, height = 5.5, bg = "white")
ggsave(paste0(OUT_STEM, ".png"), panel, width = 12, height = 5.5, dpi = 150,
       bg = "white")
stopifnot(file.exists(paste0(OUT_STEM, ".pdf")),
          file.exists(paste0(OUT_STEM, ".png")))
cat(sprintf("  Saved %s.{pdf,png}\n", OUT_STEM))

# =============================================================================
# STEP 4: Verification
# =============================================================================

cat("\n=== Verification ===\n")
cat(sprintf("chr21 median: uncorrected %.3f -> corrected %.3f (shift %.3f)\n",
            stats[chromosome == "chr21" & scale == "uncorrected", median],
            stats[chromosome == "chr21" & scale == "ploidy-corrected", median],
            shift[chromosome == "chr21", median_shift]))
cat(sprintf("%s median: uncorrected %.3f -> corrected %.3f (max |shift| %.2g)\n",
            CONTROL_CHR,
            stats[chromosome == CONTROL_CHR & scale == "uncorrected", median],
            stats[chromosome == CONTROL_CHR & scale == "ploidy-corrected", median],
            ctrl_shift))

writeLines(capture.output(sessionInfo()),
           paste0(OUT_STEM, "_session_info.txt"))
stopifnot(file.exists(paste0(OUT_STEM, "_session_info.txt")))

cat("\n=== Distribution panel complete ===\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-09-04  REBUILT around the ploidy-correction effect. The panel showed
#             only the corrected chr21 distribution against a baseMean-matched
#             non-chr21 sample, which cannot show what the correction did; it
#             now overlays uncorrected and ploidy-corrected log2FC for chr21,
#             with chr22 on both scales as the control (the correction must
#             leave it unchanged). Dropped the per-lane |robust z| scatter.
#             Input is now deseq2_all_genes_both_analyses.csv (both scales in
#             one table); the lane table is no longer read. Data preparation
#             moved to scripts/lib/ploidy_distributions.R (unit-tested).
#             Outputs renamed from chr21_vs_genome_distribution.* to
#             ploidy_correction_distributions.*; a stats CSV is written.
#
# 2026-09-01  REPLACED deviation_vs_cohort_sd with abs(dev_z) in the (since
#             removed) per-lane scatter; ADDED existence checks after ggsave.
# 2026-09-04  WIDENED the gene set to TARGET_BIOTYPES (scripts/lib/biotypes.R:
#             protein-coding, lncRNA, pseudogene; passed to
#             ploidy_distribution_long(biotypes =)) to match scripts 02 and 04.
