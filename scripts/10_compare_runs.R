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
# Usage:  Rscript scripts/10_compare_runs.R --run adjusted --against baseline
# Outputs (under the --run run):
#   tables/run_comparison_<run>_vs_<against>.csv
#   tables/lane_transitions_<run>_vs_<against>.csv
#   figures/run_comparison_<run>_vs_<against>.{pdf,png}

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("scripts/lib/run.R")
args    <- commandArgs(trailingOnly = TRUE)
run     <- load_run(args)
this    <- run$name
against <- if ("--against" %in% args) args[match("--against", args) + 1] else "baseline"
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

lab_b <- if (this == against) paste0(this, " (2)") else this
long <- rbind(dev[, .(Gene_name, run_label = against, lfc = get(lfc_a))],
              dev[, .(Gene_name, run_label = lab_b, lfc = get(lfc_b))])
long[, Gene_name := factor(Gene_name, levels = dev$Gene_name)]
long[, run_label := factor(run_label, levels = c(against, lab_b))]
cols_run <- setNames(c("#C45A00", "#762A83"), c(against, lab_b))
th <- theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 10),
                                       panel.grid.major.y = element_blank(), legend.position = "bottom")
p1 <- ggplot(long, aes(x = lfc, y = Gene_name, colour = run_label)) +
  geom_line(aes(group = Gene_name), colour = "grey75") + geom_point(size = 2, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dotted") + scale_colour_manual(values = cols_run, name = NULL) +
  labs(title = "Ploidy-corrected log2FC by run", x = "norm_log2FC", y = NULL) + th
fig <- p1
stem <- run$figure(sprintf("run_comparison_%s_vs_%s", this, against))
ggsave(paste0(stem, ".pdf"), fig, width = 6, height = 7)
ggsave(paste0(stem, ".png"), fig, width = 6, height = 7, dpi = 200)
writeLines(capture.output(sessionInfo()), paste0(stem, "_session_info.txt"))
cat(sprintf("\n  Saved: %s.{pdf,png}\n=== Comparison complete ===\n", stem))
