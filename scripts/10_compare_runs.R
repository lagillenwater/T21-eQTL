# 10_compare_runs.R
#
# Purpose: Per-gene comparison of two runs' lane tables: lane, tier,
#          ploidy-corrected log2FC, eQTL call,
#          for every gene deviating in either run, plus a lane-transition
#          table for all chr21 genes. Replaces the retired
#          10_composition_adjusted_neighborhood.R (scripts/archive/), which
#          approximated the adjustment on a log-CPM scale; here both runs are
#          full DESeq2 fits.
#
#          When both runs have the deviation-bound, base-rate and
#          matched-effect-size tables (scripts 13, 16, 18), also the
#          run-level statistics behind the negative conclusion, side by side
#          (scripts/lib/run_summary.R).
#
# Usage:  Rscript scripts/10_compare_runs.R --run adjusted --against baseline
# Outputs (under the --run run):
#   tables/run_comparison_<run>_vs_<against>.csv
#   tables/lane_transitions_<run>_vs_<against>.csv
#   tables/negative_conclusion_<run>_vs_<against>.csv   (when 13, 16, 18 ran for both)
#   figures/run_comparison_<run>_vs_<against>.{pdf,png}

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("scripts/lib/run.R"); source("scripts/lib/run_summary.R")
source("scripts/lib/figure_style.R")   # record_titles: no in-plot title
args    <- commandArgs(trailingOnly = TRUE)
run     <- load_run(args)
this    <- run$name
j       <- match("--against", args)
if (!is.na(j) && (j == length(args) || startsWith(args[j + 1], "--")))
  stop("--against needs a run name, e.g. --against baseline", call. = FALSE)
against <- if (is.na(j)) "baseline" else args[j + 1]
base    <- load_run(c("--run", against))
cat(sprintf("=== Run comparison: %s vs %s ===\n\n", this, against))

cols <- c("Gene_name", "Gene_type", "sig_lane", "tier", "norm_log2FC", "norm_padj",
          "eqtl_lane", "q_gene_bh")
a <- fread(base$table("chr21_lane_assignments.csv"))[, ..cols]
b <- fread(run$table("chr21_lane_assignments.csv"))[, ..cols]
sa <- paste0("_", against); sb <- paste0("_", this)
if (sa == sb) sb <- paste0(sb, "_2")   # self-comparison (smoke test) keeps columns distinct
m <- merge(a, b, by = c("Gene_name", "Gene_type"), suffixes = c(sa, sb), all = TRUE)
lane_a <- paste0("sig_lane", sa); lane_b <- paste0("sig_lane", sb)
lfc_a  <- paste0("norm_log2FC", sa); lfc_b <- paste0("norm_log2FC", sb)

dev <- m[get(lane_a) %in% c("DE_high", "DE_low") | get(lane_b) %in% c("DE_high", "DE_low")]
dev[, attenuation := 1 - get(lfc_b) / get(lfc_a)]
setorderv(dev, lfc_a)
fwrite(dev, run$table(sprintf("run_comparison_%s_vs_%s.csv", this, against)))

trans <- m[, .N, by = .(from = get(lane_a), to = get(lane_b))][order(from, to)]
setnames(trans, c("from", "to"), c(paste0("lane", sa), paste0("lane", sb)))
fwrite(trans, run$table(sprintf("lane_transitions_%s_vs_%s.csv", this, against)))
cat("Lane transitions (all chr21 genes):\n"); print(as.data.frame(trans))
cat("\nDeviating genes in either run:\n")
print(as.data.frame(dev[, c("Gene_name", lane_a, lane_b, lfc_a, lfc_b, "attenuation",
                            paste0("eqtl_lane", sa), paste0("eqtl_lane", sb)), with = FALSE]))

# Run names in words on the figure; an unknown run keeps its name.
RUN_WORDS <- c(baseline = "Unadjusted", adjusted = "Adjusted",
               adjusted_interaction = "Adjusted, karyotype x composition")
run_words <- function(r) if (r %in% names(RUN_WORDS)) RUN_WORDS[[r]] else r
lab_a <- run_words(against)
lab_b <- if (this == against) paste0(run_words(this), " (2)") else run_words(this)
long <- rbind(dev[, .(Gene_name, run_label = lab_a, lfc = get(lfc_a))],
              dev[, .(Gene_name, run_label = lab_b, lfc = get(lfc_b))])
long[, Gene_name := factor(Gene_name, levels = dev$Gene_name)]
long[, run_label := factor(run_label, levels = c(lab_a, lab_b))]
cols_run <- setNames(c("#C45A00", "#762A83"), c(lab_a, lab_b))
th <- theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 10),
                                       panel.grid.major.y = element_blank(), legend.position = "bottom")
p1 <- ggplot(long, aes(x = lfc, y = Gene_name, colour = run_label)) +
  geom_line(aes(group = Gene_name), colour = "grey75") + geom_point(size = 2, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dotted") + scale_colour_manual(values = cols_run, name = NULL) +
  labs(title = "Ploidy-corrected log2FC by run", x = "Ploidy-corrected log2 fold change (T21 vs D21)", y = NULL) +
  th + theme(axis.text.y = element_text(face = "italic"))
stem <- run$figure(sprintf("run_comparison_%s_vs_%s", this, against))
fig <- strip_titles(p1, basename(stem), "", file.path(run$dirs$figures, "removed_titles.txt"))
ggsave(paste0(stem, ".pdf"), fig, width = 6, height = 7)
ggsave(paste0(stem, ".png"), fig, width = 6, height = 7, dpi = 200)
writeLines(capture.output(sessionInfo()), paste0(stem, "_session_info.txt"))
outputs <- c(run$table(sprintf("run_comparison_%s_vs_%s.csv", this, against)),
             run$table(sprintf("lane_transitions_%s_vs_%s.csv", this, against)),
             paste0(stem, c(".pdf", ".png", "_session_info.txt")))

# Run-level statistics behind the negative conclusion, when both runs have them.
SUMMARY_TABLES <- c(bound_genes = "af_deviation_bound_genes.csv", background = "af_reference_background.csv",
                    base_rate = "base_rate_summary.csv", detection = "detection_vs_effect_size_summary.csv",
                    ratio = "top_variant_effect_ratio_summary.csv")
run_stats <- function(r) {
  paths <- vapply(SUMMARY_TABLES, r$table, character(1))
  if (!all(file.exists(paths))) return(NULL)
  do.call(negative_conclusion_summary, lapply(paths, fread))
}
sa_stats <- run_stats(base); sb_stats <- run_stats(run)
if (!is.null(sa_stats) && !is.null(sb_stats)) {
  neg <- merge(sa_stats, sb_stats, by = "statistic", suffixes = c(sa, sb), sort = FALSE)
  value_cols <- grep("^value_", names(neg), value = TRUE)
  neg[, (value_cols) := lapply(.SD, signif, 4), .SDcols = value_cols]
  neg_path <- run$table(sprintf("negative_conclusion_%s_vs_%s.csv", this, against))
  fwrite(neg, neg_path)
  outputs <- c(outputs, neg_path)
  cat("\nNegative-conclusion statistics by run:\n"); print(as.data.frame(neg))
} else {
  cat("\n  Skipping the negative-conclusion statistics: scripts 13, 16 and 18 have not run for both runs.\n")
}
if (!all(file.exists(outputs))) stop("failed to write: ", paste(outputs[!file.exists(outputs)], collapse = ", "))
cat(sprintf("\n  Saved: %s.{pdf,png}\n=== Comparison complete ===\n", stem))
