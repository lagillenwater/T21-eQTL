# 12_flow_map_figure.R
#
# Purpose: One two-panel figure for the manuscript: the flow of the chr21
#          genes through the dosage classification to the eQTL test (A), and
#          the deviating genes on chr21 (B). Both panels are drawn here from
#          the run's tables, so the figure always matches the run.
#
#   A  Every classified chr21 gene (the lane table) through three stages:
#      expected dosage / outside the expectation / not assessable; then the
#      sub-lanes (not DE, higher than expected, lower than expected; low
#      expression, high repeats); then, for the deviating genes, the
#      gene-level cis-eQTL test (detected, tested and not detected, no
#      variant to test). Flows end where a gene's path ends. Drawn with
#      ggsankey. Colours follow the figure palette: red and blue for the
#      direction of deviation, filled / open / grey for the eQTL outcome,
#      greys for the dosage lanes.
#   B  The chr21 map of the deviating genes, deviating_map_plot() in
#      scripts/lib/eqtl_figures.R (the same map script 11 writes alone).
#
# Inputs (per run): tables/chr21_lane_assignments.csv (script 04),
#   processed/eqtl_supported_genes.csv (script 02), data/chr21_gene_positions.csv,
#   data/Whole_Blood.v10.eGenes.txt.gz
# Outputs (per run): figures/chr21_flow_and_map.{pdf,png}
#
# Usage: Rscript scripts/12_flow_map_figure.R --run adjusted
# Requires ggsankey (GitHub: davidsjoberg/ggsankey; install_packages.R).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggsankey)
  library(patchwork)
})
source("scripts/lib/eqtl_figures.R")   # deviating_map_data, deviating_map_plot
source("scripts/lib/figure_style.R")   # record_titles
source("scripts/lib/run.R"); run <- load_run()
TITLES <- file.path(run$dirs$figures, "removed_titles.txt")
EGENES <- "data/Whole_Blood.v10.eGenes.txt.gz"
COL_DIR   <- c(up = "#B2182B", down = "#2166AC")
CHR21_LEN <- 46709983                    # GRCh38
CEN       <- c(10864561, 12915808)

cat(sprintf("=== T21-eQTL: lane flow + chr21 map [run: %s] ===\n\n", run$name))
lanes <- fread(run$table("chr21_lane_assignments.csv"))

# =============================================================================
# A: lane flow
# =============================================================================
# One row per gene, one column per stage; NA where the gene's path ends.
ROOT <- "HSA21 genes"
NODE_ORDER <- c(ROOT,
                "Expected dosage", "Outside the expectation", "Not assessable",
                "Did not reach significance", "Higher than expected", "Lower than expected",
                "Low expression", "High repeats", "Not estimable",
                "cis-eQTL detected", "tested, not detected", "no variant to test")
NODE_FILL <- c("grey40",
               "grey65", "grey45", "grey80",
               "grey60", COL_DIR[["up"]], COL_DIR[["down"]],
               "grey85", "grey85", "grey85",
               "grey15", "white", "grey70")
names(NODE_FILL) <- NODE_ORDER

g <- lanes[, .(sig_lane, eqtl_lane)]
g[, `:=`(
  s1 = ROOT,
  s2 = fcase(sig_lane == "Expected_dosage", "Expected dosage",
             sig_lane %in% c("DE_high", "DE_low", "Not_DE_outside_noise"), "Outside the expectation",
             sig_lane %in% c("Low_expression", "High_repeats", "Not_assessable"), "Not assessable"),
  s3 = fcase(sig_lane == "Not_DE_outside_noise", "Did not reach significance",
             sig_lane == "DE_high", "Higher than expected",
             sig_lane == "DE_low", "Lower than expected",
             sig_lane == "Low_expression", "Low expression",
             sig_lane == "High_repeats", "High repeats",
             sig_lane == "Not_assessable", "Not estimable"),
  s4 = fcase(sig_lane %in% c("DE_high", "DE_low") & eqtl_lane == "cis_eqtl", "cis-eQTL detected",
             sig_lane %in% c("DE_high", "DE_low") & eqtl_lane == "no_cis_eqtl", "tested, not detected",
             sig_lane %in% c("DE_high", "DE_low") & eqtl_lane == "no_GTEx_data", "no variant to test"))]
if (anyNA(g$s2)) stop("unknown sig_lane: ", paste(unique(g[is.na(s2), sig_lane]), collapse = ", "))
if (g[sig_lane %in% c("DE_high", "DE_low"), anyNA(s4)]) stop("deviating gene without an eQTL outcome")

long <- as.data.table(make_long(as.data.frame(g), s1, s2, s3, s4))
long <- long[!is.na(node)]
counts <- long[, .N, by = .(x, node)]
cat("  Flow nodes:\n"); print(counts[, .(stage = x, node, genes = N)])
long <- merge(long, counts, by = c("x", "node"))
long[, `:=`(node = factor(node, levels = rev(NODE_ORDER)),
            next_node = factor(next_node, levels = rev(NODE_ORDER)),
            label = sprintf("%s  %d", node, N))]

pA <- ggplot(long, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                       fill = node, label = label)) +
  # space: the gap between nodes, in genes. Wide enough that the labels of the
  # small nodes (6, 7, 11 genes; 2, 3, 8) clear each other at a page width.
  geom_sankey(flow.alpha = 0.45, node.color = "grey30", width = 0.08, space = 24) +
  # Labels on a 55% white box: several sit where flows leave their node, so
  # the text stays readable and the flows behind it still show.
  geom_sankey_label(size = 2.7, hjust = 0, colour = "grey10", fill = grDevices::adjustcolor("white", alpha.f = 0.55),
                    border.colour = NA, space = 24, position = position_nudge(x = 0.06)) +
  scale_fill_manual(values = NODE_FILL, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.1, 1.0))) +
  theme_void(base_size = 9)

# =============================================================================
# B: chr21 map of the deviating genes
# =============================================================================
gtex_tested_ids <- sub("[.].*$", "", fread(EGENES, select = "gene_id")$gene_id)
bands <- deviating_map_data(lanes, fread(run$processed("eqtl_supported_genes.csv")),
                            fread("data/chr21_gene_positions.csv"), gtex_tested_ids)
pB <- deviating_map_plot(bands, COL_DIR, CHR21_LEN, CEN)
record_titles("chr21_flow_and_map", "B", "Deviating chr21 genes and their cis-eQTL outcome",
              sprintf(paste0("Bands at the TSS of each of the %d deviating genes (GRCh38); labels above the bar are ",
                             "higher, below are lower; overlapping bands share the bar height. Centromere in dark grey."),
                      nrow(bands)), TITLES)

fig <- pA / pB + plot_layout(heights = c(1.3, 0.8)) + plot_annotation(tag_levels = "A")
# The two panels again, for the combined Figure 2 (script 20).
save_panels(list(flow = pA, map = pB), "chr21_flow_and_map", run)
stem <- run$figure("chr21_flow_and_map")
ggsave(paste0(stem, ".pdf"), fig, width = 10, height = 9.5)
ggsave(paste0(stem, ".png"), fig, width = 10, height = 9.5, dpi = 200)
for (f in paste0(stem, c(".pdf", ".png"))) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat(sprintf("  Saved: %s\n", f))
}
cat("\n=== Combined figure complete ===\n")
