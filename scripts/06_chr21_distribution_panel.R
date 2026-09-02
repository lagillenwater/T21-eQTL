# 15_chr21_distribution_panel.R
#
# Purpose: Reproduce the paper's Fig 2/3-style chromosome-21 vs genome
#          distribution comparison for our cohort. This is the headline
#          population-level evidence that should sit BEFORE any per-gene
#          drill-down: it shows whether chromosome 21's ploidy-normalized
#          fold-change distribution differs from baseMean-matched
#          non-chr21 cohort variation. Scope what we can claim:
#            - If chr21 distribution overlaps non-chr21, the paper's
#              "no dosage compensation" finding holds at population level
#              for our cohort.
#            - It does NOT rule out per-gene compensation; the per-gene
#              drill-down (script 04) is still needed.
#
# Inputs:
#   - results/tables/deseq2_all_genes_ploidy_normalized.csv
#   - results/tables/chr21_lane_assignments.csv  (for annotation overlay)
#
# Outputs:
#   - results/figures/chr21_vs_genome_distribution.pdf
#   - results/figures/chr21_vs_genome_distribution.png
#   - results/figures/chr21_vs_genome_distribution_session_info.txt
#
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# Match the rest of the pipeline (scripts 09 and 12) - apples-to-apples
# comparison requires both sides of the chr21-vs-genome distribution to
# be on the same gene-type subset.
RESTRICT_TO_PROTEIN_CODING <- TRUE

cat("=== T21-eQTL: Chr21 vs Genome Distribution Panel ===\n\n")

# =============================================================================
# STEP 1: Load and partition log2FC distributions
# =============================================================================

cat("Step 1: Loading ploidy-normalized DESeq2 results...\n")

all_lfc <- fread("results/tables/deseq2_all_genes_ploidy_normalized.csv")
all_lfc <- all_lfc[!is.na(log2FoldChange) & !is.na(baseMean)]
cat(sprintf("  Genes with non-NA log2FC: %d\n", nrow(all_lfc)))

if (RESTRICT_TO_PROTEIN_CODING) {
  n_before <- nrow(all_lfc)
  all_lfc <- all_lfc[Gene_type == "protein_coding"]
  cat(sprintf("  Restricted to protein-coding: %d -> %d\n",
              n_before, nrow(all_lfc)))
}

chr21 <- all_lfc[Chr == "chr21"]
other <- all_lfc[Chr != "chr21"]

cat(sprintf("  chr21: %d   non-chr21: %d\n", nrow(chr21), nrow(other)))

# baseMean-matched control set: trim non-chr21 to chr21's baseMean range
# AND match the chr21 baseMean distribution by quantile binning. Keeps the
# comparison fair against expression-level confounding.
bm_q <- quantile(chr21$baseMean, c(0.05, 0.95))
ctrl_pool <- other[baseMean >= bm_q[1] & baseMean <= bm_q[2]]
cat(sprintf("  Control pool (non-chr21, baseMean-trimmed): %d\n",
            nrow(ctrl_pool)))

# Quantile-stratified random sample so the control's baseMean distribution
# matches chr21's at the same n.
chr21[, bm_decile := cut(baseMean,
                         breaks = quantile(baseMean,
                                           probs = seq(0, 1, 0.1),
                                           na.rm = TRUE),
                         include.lowest = TRUE,
                         labels = FALSE)]
ctrl_pool[, bm_decile := cut(baseMean,
                             breaks = quantile(chr21$baseMean,
                                               probs = seq(0, 1, 0.1),
                                               na.rm = TRUE),
                             include.lowest = TRUE,
                             labels = FALSE)]
target_per_decile <- chr21[, .N, by = bm_decile]
matched <- ctrl_pool[!is.na(bm_decile)][
  , .SD[sample(.N, min(.N, target_per_decile[bm_decile == .BY[[1]], N]))],
  by = bm_decile]
cat(sprintf("  baseMean-matched control sample: %d\n", nrow(matched)))

# =============================================================================
# STEP 2: Summary statistics
# =============================================================================

stat_block <- function(x, label) {
  data.table(
    set    = label,
    n      = length(x),
    median = median(x),
    mean   = mean(x),
    sd     = sd(x),
    p05    = quantile(x, 0.05),
    p95    = quantile(x, 0.95)
  )
}
stats <- rbindlist(list(
  stat_block(chr21$log2FoldChange,   "chr21 (n=322)"),
  stat_block(matched$log2FoldChange, "non-chr21 (matched)"),
  stat_block(other$log2FoldChange,   "non-chr21 (all)")
))
cat("\n=== Distribution stats ===\n")
print(stats)

ks <- ks.test(chr21$log2FoldChange, matched$log2FoldChange)
cat(sprintf("\nKS test chr21 vs matched non-chr21: D=%.3f  p=%.3g\n",
            ks$statistic, ks$p.value))

# =============================================================================
# STEP 3: Density + ECDF panels
# =============================================================================

cat("\nStep 3: Building plots...\n")

dens_df <- rbindlist(list(
  data.table(set = "chr21",            log2FC = chr21$log2FoldChange),
  data.table(set = "non-chr21 matched", log2FC = matched$log2FoldChange)
))

p_density <- ggplot(dens_df, aes(x = log2FC, fill = set, color = set)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values  = c("chr21" = "#D62728",
                                "non-chr21 matched" = "#1F77B4")) +
  scale_color_manual(values = c("chr21" = "#D62728",
                                "non-chr21 matched" = "#1F77B4")) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(
    title    = "Density of ploidy-normalized log2 fold change",
    subtitle = sprintf(paste0("chr21 (n=%d) vs baseMean-matched non-chr21 ",
                              "(n=%d). KS p=%.3g"),
                       nrow(chr21), nrow(matched), ks$p.value),
    x = "log2(T21 / Control)  [chr21 corrected for ploidy]",
    y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"))

p_ecdf <- ggplot(dens_df, aes(x = log2FC, color = set)) +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("chr21" = "#D62728",
                                "non-chr21 matched" = "#1F77B4")) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(
    title = "ECDF: chr21 vs non-chr21 (matched)",
    x     = "log2(T21 / Control)",
    y     = "Cumulative fraction",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# =============================================================================
# STEP 4: Per-lane magnitude annotation panel
# =============================================================================

lanes <- fread("results/tables/chr21_lane_assignments.csv")
lanes_de <- lanes[sig_lane %in% c("DE_low", "DE_high") &
                  low_expr == FALSE & high_repeat == FALSE]
lanes_de[, lane_label := fcase(
  eqtl_lane == "cis_eqtl",     "cis-eQTL detected",
  eqtl_lane == "no_cis_eqtl",  "eQTL-tested, none detected",
  eqtl_lane == "no_GTEx_data", "no GTEx eQTL data")]
lanes_de[, lane_label := factor(lane_label,
  levels = c("cis-eQTL detected",
             "eQTL-tested, none detected",
             "no GTEx eQTL data"))]

ctrl_sd <- sd(matched$log2FoldChange)

p_lane <- ggplot(lanes_de,
                 aes(x = lane_label, y = abs(dev_z),
                     color = sig_lane)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = c(1, 2), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("DE_low" = "#FF7F0E",
                                "DE_high" = "#1F77B4")) +
  labs(
    title    = "Per-gene deviation magnitude vs cohort noise, by lane",
    subtitle = paste0("Dashed lines at |robust z| = 1 and 2 (chr21 median/MAD null). ",
                      "Genes <1 SD are within typical cohort variation; ",
                      "|z| > 2 sit clearly outside it."),
    x        = NULL,
    y        = "|robust z| (chr21 median/MAD units)",
    color    = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),
        legend.position = "top")

# =============================================================================
# STEP 5: Compose and save
# =============================================================================

cat("\nStep 5: Composing panel...\n")

panel <- (p_density | p_ecdf) / p_lane +
  plot_annotation(
    title = paste0("chr21 ploidy-normalized fold-change distribution and ",
                   "per-gene magnitude triage"),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("results/figures/chr21_vs_genome_distribution.pdf", panel,
       width = 13, height = 9, bg = "white")
ggsave("results/figures/chr21_vs_genome_distribution.png", panel,
       width = 13, height = 9, dpi = 150, bg = "white")
stopifnot(file.exists("results/figures/chr21_vs_genome_distribution.pdf"),
          file.exists("results/figures/chr21_vs_genome_distribution.png"))
cat("  Saved chr21_vs_genome_distribution.{pdf,png}\n")

# =============================================================================
# STEP 6: Verification
# =============================================================================

cat("\n=== Verification ===\n")
cat(sprintf("chr21 SD: %.3f   non-chr21 matched SD: %.3f   ratio: %.2f\n",
            sd(chr21$log2FoldChange), ctrl_sd,
            sd(chr21$log2FoldChange) / ctrl_sd))
cat(sprintf("chr21 median: %.3f   non-chr21 matched median: %.3f   shift: %.3f\n",
            median(chr21$log2FoldChange), median(matched$log2FoldChange),
            median(chr21$log2FoldChange) - median(matched$log2FoldChange)))

cat("\nLane-magnitude triage (n with |robust z| >= 2):\n")
print(lanes_de[, .(n_total = .N,
                   n_above_1z = sum(abs(dev_z) > 1),
                   n_above_2z = sum(abs(dev_z) > 2)),
               by = .(sig_lane, lane_label)])

writeLines(capture.output(sessionInfo()),
           "results/figures/chr21_vs_genome_distribution_session_info.txt")
stopifnot(file.exists("results/figures/chr21_vs_genome_distribution_session_info.txt"))

cat("\n=== Distribution panel complete ===\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-09-01  REPLACED deviation_vs_cohort_sd (dropped from the lane table by
#             the 2026-08-31 threshold change) with abs(dev_z), the robust z
#             against the chr21 median/MAD null. abs() is required: the old
#             column was a magnitude, dev_z is signed, and a naive rename made
#             the DE_low tallies read 0. Labels updated from "cohort-noise SDs"
#             to robust-z units. This script was omitted from the branch's
#             full-chain check, which is how the break survived review.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
#
# 2026-09-01  ADDED existence checks for both figure files directly after the
#             ggsave() calls; previously only the session-info file was
#             verified (CodeRabbit review, PR #3).
