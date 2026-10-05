# S3_distal_inflation.R
#
# Purpose: Supplement figure for genotype structure in the within-T21 eQTL
#          tests: the distal-variant inflation scan (script 17) of two runs
#          side by side, and the run's covariates against chr21 genotype.
#
#   Variants on the long arm only (script 17's "q-arm only" set, >= 13 Mb):
#   the short arm and pericentromere carry GRCh38's false duplications of
#   21q22.3, where a variant's coordinate puts a q-arm gene's own cis signal
#   more than 5 Mb away (GATD3A, PWP2, SIK1, ICOSLG reach p < 1e-40 that
#   way), so the distance rule cannot be applied there. Script 17 reports
#   the all-variant numbers beside these.
#
#   A. QQ of the distal tests, HSA21 assessable genes and the HSA22 reference
#      genes, each run with its subject-permutation 95% envelope.
#   B. lambda per gene set with its permutation 95% range, both runs.
#   C. The --run run's covariates, each regressed on every common chr21
#      variant: a covariate that tracks genotype across the chromosome
#      carries the structure, and adjusting for it removes the structure by
#      proxy.
#
#   Descriptive: lambda and the QQ curve say whether the p-values depart from
#   uniform, not what causes the departure.
#
# Usage:   Rscript scripts/S3_distal_inflation.R --run adjusted --against baseline
# Inputs:  tables/distal_inflation_{qq,summary,covariates}.csv of both runs
#          (script 17; covariates only for a run with covariates)
# Outputs: results/runs/<run>/figures/S3_distal_inflation.{pdf,png}
#          results/runs/<run>/tables/S3_distal_inflation_lambda.csv

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("scripts/lib/covariates.R")    # display_covariate
source("scripts/lib/run.R")
args    <- commandArgs(trailingOnly = TRUE)
run     <- load_run(args)
j       <- match("--against", args)
if (!is.na(j) && (j == length(args) || startsWith(args[j + 1], "--")))
  stop("--against needs a run name, e.g. --against baseline", call. = FALSE)
against <- if (is.na(j)) "baseline" else args[j + 1]
base    <- load_run(c("--run", against))
cat(sprintf("=== S3 distal-variant inflation: %s vs %s ===\n\n", run$name, against))

RUN_LAB  <- setNames(c("Unadjusted", "Adjusted"), c(against, run$name))
VARIANTS <- "q-arm only"                                   # long-arm variants (>= 13 Mb)
COL_RUN  <- setNames(c("#C45A00", "#762A83"), RUN_LAB)   # as in 10_compare_runs.R
SET_LAB  <- c("chr21 assessable" = "HSA21 assessable genes", "chr21 deviating" = "HSA21 deviating genes",
              "chr22 reference" = "HSA22 reference genes")
th <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), panel.grid.minor = element_blank(),
        legend.position = "bottom")

read_both <- function(file) rbind(
  fread(base$table(file))[, run_label := RUN_LAB[[against]]],
  fread(run$table(file))[, run_label := RUN_LAB[[run$name]]])

# ---- A: QQ per gene set, stacked at full width, zoomed, full-range inset ---
# Inflation of a few percent lifts the whole curve by a few hundredths of a
# -log10 p unit, invisible over the full range, so each panel shows expected
# 0 to ZOOM_X, where the unadjusted curve leaves its permutation range; the
# full range sits in an inset with the zoomed window outlined.
ZOOM_X <- 4
qq <- read_both("distal_inflation_qq.csv")[variant_set == VARIANTS & gene_set %in% names(SET_LAB)]
qq[, `:=`(set_f = factor(SET_LAB[gene_set], levels = SET_LAB),
          run_f = factor(run_label, levels = RUN_LAB))]
X_MAX <- max(qq$expected) * 1.02
Y_MAX <- max(qq$perm_hi, qq$observed) * 1.02
qq_layers <- function(d) list(
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4),
  geom_ribbon(data = d, aes(ymin = perm_lo, ymax = perm_hi, fill = run_f), colour = NA, alpha = 0.15),
  geom_line(data = d, linewidth = 0.6),
  scale_colour_manual(values = COL_RUN, name = NULL),
  scale_fill_manual(values = COL_RUN, guide = "none"))
qq_panel <- function(s, title, subtitle = NULL) {
  d <- qq[set_f == s]
  zoom_y <- max(d[expected <= ZOOM_X, c(observed, perm_hi)]) * 1.02
  main <- ggplot(d[expected <= ZOOM_X * 1.05], aes(expected, observed, colour = run_f)) +
    qq_layers(d[expected <= ZOOM_X * 1.05]) +
    coord_cartesian(xlim = c(0, ZOOM_X), ylim = c(0, zoom_y), expand = FALSE) +
    labs(x = expression(Expected ~ -log[10] ~ italic(p)), y = expression(Observed ~ -log[10] ~ italic(p)),
         title = title, subtitle = subtitle) +
    th
  inset <- ggplot(d, aes(expected, observed, colour = run_f)) +
    qq_layers(d) +
    annotate("rect", xmin = 0, xmax = ZOOM_X, ymin = 0, ymax = zoom_y, fill = NA, colour = "grey30",
             linewidth = 0.3) +
    coord_cartesian(xlim = c(0, X_MAX), ylim = c(0, Y_MAX), expand = FALSE) +
    guides(colour = "none") +
    labs(x = NULL, y = NULL, title = "full range") +
    theme_bw(base_size = 6) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 6),
          plot.background = element_rect(colour = "grey70"))
  main + inset_element(inset, left = 0.70, bottom = 0.04, right = 0.99, top = 0.62)
}
pA <- qq_panel(SET_LAB[["chr21 assessable"]],
               "A  HSA21 assessable genes",
               "T21 expression against common HSA21 long-arm variants >= 5 Mb away; shaded: permutation 95% range") /
  qq_panel(SET_LAB[["chr22 reference"]], "HSA22 reference genes")

# ---- B: lambda per gene set --------------------------------------------------
lam <- read_both("distal_inflation_summary.csv")[variant_set == VARIANTS & gene_set %in% names(SET_LAB)]
lam[, `:=`(set_f = factor(SET_LAB[gene_set], levels = rev(SET_LAB)),
           run_f = factor(run_label, levels = RUN_LAB))]
fwrite(lam[, .(run = run_label, gene_set, n_genes, n_tests, lambda, perm_lo, perm_hi, perm_p,
               lambda_inverse_normal)], run$table("S3_distal_inflation_lambda.csv"))
# Explicit offsets rather than position_dodge, which does not move text on a
# discrete y axis: unadjusted below, adjusted above each set's row.
lam[, y := as.numeric(set_f) + fifelse(run_f == RUN_LAB[[run$name]], 0.15, -0.15)]
pB <- ggplot(lam, aes(lambda, y, colour = run_f)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(xmin = perm_lo, xmax = perm_hi), linewidth = 2.2, colour = "grey80") +
  geom_point(size = 2) +
  geom_text(aes(x = pmax(lambda, perm_hi), label = sprintf("n = %d", n_genes)), hjust = -0.3,
            size = 2.4, colour = "grey30", show.legend = FALSE) +
  scale_y_continuous(breaks = seq_along(levels(lam$set_f)), labels = levels(lam$set_f),
                     expand = expansion(add = 0.4)) +
  scale_colour_manual(values = COL_RUN, guide = "none") +   # the line legend of A and B names the runs
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = expression(lambda ~ "(genomic inflation)"), y = NULL,
       title = "B  Inflation by gene set",
       subtitle = "Grey: permutation 95% range (2.5th-97.5th percentile)") +
  th + theme(panel.grid.major.y = element_blank())

# pA's two stacked panels flatten into the outer layout, so heights are per row.
panels <- pA / pB
heights <- c(0.8, 0.8, 0.6)

# ---- C: the run's covariates against chr21 genotype --------------------------
cov_path <- run$table("distal_inflation_covariates.csv")
if (file.exists(cov_path)) {
  cv <- fread(cov_path)
  cv[, cov_f := factor(display_covariate(covariate), levels = rev(display_covariate(covariate)))]
  pC <- ggplot(cv, aes(lambda, cov_f)) +
    geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(xmin = perm_lo, xmax = perm_hi), linewidth = 2.2, colour = "grey80") +
    # neutral ink: the run colours of A-C would read as a run comparison here
    geom_point(aes(shape = above_range), size = 1.8, colour = "grey15") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(x = expression(lambda ~ "(covariate on each common HSA21 variant)"), y = NULL,
         title = "C  Covariates against HSA21 genotype",
         subtitle = "Adjusted run; filled: outside its permutation 95% range") +
    th
  panels <- pA / (pB | pC)
  heights <- c(0.8, 0.8, 1)
}

fig <- panels + plot_layout(heights = heights, guides = "collect") & theme(legend.position = "bottom")
stem <- run$figure("S3_distal_inflation")
ggsave(paste0(stem, ".pdf"), fig, width = 9, height = 13)
ggsave(paste0(stem, ".png"), fig, width = 9, height = 13, dpi = 300)
print(lam[, .(run_label, gene_set, n_genes, lambda, perm_lo, perm_hi, perm_p)])
cat("\nDone.\n")
