# 11_eqtl_figures.R
#
# Purpose: Three figures for the within-T21 cis-eQTL stage, drawn from tables
#          scripts 03 and 04 already wrote. No statistic is recomputed here.
#
#   eqtl_dosage_panels   Expression by alt-allele dosage in T21, box + jitter,
#                        one panel per gene, in four groups: DE high, DE low,
#                        positive controls (strong GTEx eGenes on their own
#                        best variant) and negative controls (each deviating
#                        gene against the best variant of its decoy set).
#                        The variant is the best variant of the gene-level
#                        permutation test, the one the classification rests
#                        on, so detected and not-detected genes are shown alike.
#   eqtl_effect_sizes    Scatter: within-T21 slope at the best variant (x)
#                        against -log10 gene-level permutation q (y), for the
#                        deviating genes (by direction), their decoy variant
#                        sets and the positive controls, with the q = 0.05
#                        line. Also writes tables/eqtl_best_variant_effects.csv.
#   chr21_deviating_map  Schematic chr21 (GRCh38) with a band at every
#                        deviating gene's TSS, coloured by eQTL outcome; gene
#                        labels coloured by direction (red up, blue down).
#
# Inputs (per run):
#   - tables/chr21_lane_assignments.csv, eqtl_gene_level_perm.csv,
#     eqtl_control_negative.csv, eqtl_control_positive.csv
#   - processed/eqtl_supported_genes.csv (TSS), genotypes_filtered.csv,
#     expression_adjusted.csv, analysis_cohort.csv
#   - data/chr21_gene_positions.csv (TSS for deviating genes GTEx lacks)
# Outputs (per run): figures/eqtl_dosage_panels.{pdf,png},
#   figures/eqtl_effect_sizes.{pdf,png}, figures/chr21_deviating_map.{pdf,png},
#   tables/eqtl_best_variant_effects.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("scripts/lib/cohort.R")        # subject_id_from_labid
source("scripts/lib/eqtl_fit.R")      # expr_from_matrix
source("scripts/lib/eqtl_figures.R")  # panel_variants, overview_rows, map_bands
source("scripts/lib/run.R"); run <- load_run()

FDR_GENE <- run$thresholds$fdr_gene

# Direction colours match the chr21 volcano (script 07); outcome colours are
# shared by the overview and the map. Both sets pass the colour-vision
# separation checks in light mode.
COL_DIR   <- c(up = "#B2182B", down = "#2166AC")
COL_EQTL  <- c(cis_eqtl = "#1B7837", no_cis_eqtl = "#E08214", no_GTEx_data = "grey55")
LAB_EQTL  <- c(cis_eqtl = "cis-eQTL detected", no_cis_eqtl = "Tested, not detected",
               no_GTEx_data = "No GTEx coverage")
COL_GROUP <- c("DE high" = COL_DIR[["up"]], "DE low" = COL_DIR[["down"]],
               "Positive control (GTEx eGenes)" = "#1B7837",
               "Negative control (decoy variants)" = "grey45")

# GRCh38 chr21: length and centromere (UCSC cytoBand acen bands).
CHR21_LEN <- 46709983
CEN       <- c(10864561, 12915808)

cat(sprintf("=== T21-eQTL: eQTL figures [run: %s] ===\n\n", run$name))

# =============================================================================
# STEP 1: Load tables
# =============================================================================

lanes  <- fread(run$table("chr21_lane_assignments.csv"))
perm   <- fread(run$table("eqtl_gene_level_perm.csv"))
neg    <- fread(run$table("eqtl_control_negative.csv"))
pos    <- fread(run$table("eqtl_control_positive.csv"))
roster <- fread(run$processed("eqtl_supported_genes.csv"))
extra  <- fread("data/chr21_gene_positions.csv")
cat(sprintf("  deviating genes: %d; tested: %d; decoy tests: %d; positive controls: %d\n",
            sum(lanes$sig_lane %in% c("DE_high", "DE_low")),
            sum(!is.na(perm$p_gene_perm)), sum(!is.na(neg$p_gene_perm)), nrow(pos)))

save_fig <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 200, limitsize = FALSE)
  for (f in paste0(stem, c(".pdf", ".png"))) {
    if (!file.exists(f)) stop("failed to write ", f)
    cat(sprintf("  Saved: %s\n", f))
  }
}

# =============================================================================
# STEP 2: Dosage panels
# =============================================================================

cat("\nStep 2: Dosage panels...\n")
panels <- panel_variants(perm, lanes, neg, pos)
cat(sprintf("  panels: %s\n",
            paste(sprintf("%s %d", levels(panels$panel_group),
                          tabulate(panels$panel_group, nbins = nlevels(panels$panel_group))),
                  collapse = "; ")))

geno     <- fread(run$processed("genotypes_filtered.csv"))
cohort   <- run$cohort()
E        <- run$expression()
meta_t21 <- cohort[Karyotype == "T21"]
meta_t21[, subject_id := subject_id_from_labid(LabID)]

geno_t21 <- geno[karyotype == "T21" & !is.na(alt_dosage) &
                   variant_id %in% panels$variant_id &
                   subject_id %in% meta_t21$subject_id,
                 .(alt_dosage = mean(alt_dosage)), by = .(variant_id, subject_id)]
subj <- sort(unique(geno_t21$subject_id))
expr_tbl <- rbindlist(lapply(unique(panels$Gene_name), function(g)
  data.table(Gene_name = g, subject_id = subj,
             expr = expr_from_matrix(g, subj, E, meta_t21))))

plot_df <- merge(panels, geno_t21, by = "variant_id", allow.cartesian = TRUE)
plot_df <- merge(plot_df, expr_tbl, by = c("Gene_name", "subject_id"))
plot_df <- plot_df[!is.na(expr)]
panels[, facet := sprintf("%s\n%s\nq = %.2g%s", Gene_name, sub("_b38$", "", variant_id),
                          q_gene_bh, fifelse(is.na(decoy_gene), "",
                                             paste0("  (decoy: ", decoy_gene, ")")))]
plot_df <- merge(plot_df, panels[, .(Gene_name, panel_group, variant_id, facet)],
                 by = c("Gene_name", "panel_group", "variant_id"))
plot_df[, facet := factor(facet, levels = panels$facet)]
plot_df[, dosage := factor(alt_dosage, levels = 0:3)]

NCOL_PANELS <- 5
group_plot <- function(df, title, colour) {
  ng <- uniqueN(df$facet)
  ggplot(df, aes(x = dosage, y = expr)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92", colour = "grey40", linewidth = 0.35) +
    geom_jitter(width = 0.18, size = 0.55, alpha = 0.45, colour = colour) +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(~ facet, scales = "free_y", ncol = min(NCOL_PANELS, ng)) +
    labs(title = title, x = "Alt-allele dosage in T21 (0-3)",
         y = "Expression (log2-CPM, run artifact)") +
    theme_bw(base_size = 8.5) +
    theme(strip.text = element_text(size = 6.3, lineheight = 0.9),
          plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank())
}
group_plots <- lapply(levels(panels$panel_group), function(g) {
  d <- plot_df[panel_group == g]
  if (nrow(d) == 0) return(NULL)
  group_plot(d, g, COL_GROUP[[g]])
})
keep <- !vapply(group_plots, is.null, logical(1))
group_rows <- vapply(levels(panels$panel_group)[keep], function(g)
  ceiling(uniqueN(plot_df[panel_group == g, facet]) / NCOL_PANELS), numeric(1))
fig_panels <- wrap_plots(group_plots[keep], ncol = 1, heights = group_rows)
save_fig(fig_panels, run$figure("eqtl_dosage_panels"),
         width = 11, height = max(6, 2.3 * sum(group_rows) + 0.8 * sum(keep)))

# =============================================================================
# STEP 3: Per-gene effect sizes
# =============================================================================
# Within-T21 slope of expression on alt-allele dosage at each gene's best
# variant (own cis variants), at the best variant of its decoy set (negative
# control), and for the positive controls. The control tables store only p,
# so the single best variant is refit here with fit_variants(), the fit script
# 03 used; the refit p must reproduce the stored min_p_obs, which ties every
# effect size to the test it came from. The alt allele is arbitrary, so the
# figure shows the absolute slope with its 95% CI.

cat("\nStep 3: Effect sizes...\n")
effects <- best_variant_effects(plot_df)

stored <- rbindlist(list(
  perm[!is.na(best_variant), .(Gene_name, panel_group = NA_character_, variant_id = best_variant, min_p_obs)],
  pos[!is.na(best_variant),  .(Gene_name, panel_group = PANEL_GROUPS[3], variant_id = best_variant, min_p_obs)],
  neg[!is.na(best_variant),  .(Gene_name, panel_group = PANEL_GROUPS[4], variant_id = best_variant, min_p_obs)]))
chk <- merge(effects, stored[is.na(panel_group), .(Gene_name, variant_id, min_p_obs)],
             by = c("Gene_name", "variant_id"))[panel_group %in% PANEL_GROUPS[1:2]]
chk <- rbind(chk, merge(effects, stored[!is.na(panel_group)],
                        by = c("Gene_name", "panel_group", "variant_id")))
stopifnot(nrow(chk) == nrow(effects),
          max(abs(log10(chk$p) - log10(chk$min_p_obs))) < 1e-6)
cat(sprintf("  Refit p matches the stored best-variant p for all %d panels\n", nrow(chk)))
fwrite(effects, run$table("eqtl_best_variant_effects.csv"))
cat("  Wrote", run$table("eqtl_best_variant_effects.csv"), "\n")

es <- attach_effects(overview_rows(perm, lanes, neg, pos), effects)
es[, `:=`(own_abs = abs(own_slope), decoy_abs = abs(decoy_slope))]
es[, `:=`(own_lo = pmax(0, own_abs - 1.96 * own_se), own_hi = own_abs + 1.96 * own_se,
          decoy_lo = pmax(0, decoy_abs - 1.96 * decoy_se), decoy_hi = decoy_abs + 1.96 * decoy_se)]
print(es[, .(Gene_name, block, detected_own, own = round(own_slope, 3), own_se = round(own_se, 3),
             decoy = round(decoy_slope, 3), decoy_se = round(decoy_se, 3))])

# One scatter: x = within-T21 slope at the best variant, y = -log10 of the
# gene-level permutation q that decides detection, so the q = 0.05 line splits
# detected from not detected. Deviating genes (by direction), their decoy sets
# and the positive controls share the axes. With 1000 permutations the q has a
# floor near 1/1001, so the strongest genes share the top row and spread only
# along x. Genes with no GTEx variants have neither value; the subtitle names them.
LAB_SET <- c(up = "Higher than expected (DE_high)", down = "Lower than expected (DE_low)",
             pos = "Positive control (strong GTEx eGene)", decoy = "Decoy variant set (negative control)")
COL_SET <- setNames(c(COL_DIR[["up"]], COL_DIR[["down"]], "#1B7837", "grey55"), LAB_SET)
SHP_SET <- setNames(c(16, 16, 17, 18), LAB_SET)
untested <- sort(es[is.na(own_slope) & block != "Positive control", Gene_name])
sc <- rbindlist(list(
  es[!is.na(own_slope), .(Gene_name, slope = own_slope, q = q_own,
                          set = fcase(block == "DE high", LAB_SET[["up"]],
                                      block == "DE low",  LAB_SET[["down"]],
                                      default =           LAB_SET[["pos"]]))],
  es[!is.na(decoy_slope), .(Gene_name, slope = decoy_slope, q = q_decoy, set = LAB_SET[["decoy"]])]))
sc[, set := factor(set, levels = LAB_SET)]
sc[, y := -log10(q)]
x_lim <- max(abs(sc$slope), na.rm = TRUE) * 1.08
lab <- sc[set != LAB_SET[["decoy"]]]

fig_effects <- ggplot(sc, aes(x = slope, y = y)) +
  geom_hline(yintercept = -log10(FDR_GENE), linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(data = sc[set == LAB_SET[["decoy"]]], aes(colour = set, shape = set), size = 2.6, alpha = 0.85) +
  geom_point(data = lab, aes(colour = set, shape = set), size = 2.5) +
  geom_text_repel(data = lab, aes(label = Gene_name, colour = set), size = 2.4,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
                  segment.alpha = 0.6, box.padding = 0.35, force = 4, seed = 1,
                  show.legend = FALSE) +
  annotate("text", x = -x_lim, y = -log10(FDR_GENE), label = sprintf("q = %.2f", FDR_GENE),
           hjust = 0, vjust = -0.5, size = 2.5, colour = "grey35") +
  scale_colour_manual(values = COL_SET, drop = FALSE, name = NULL) +
  scale_shape_manual(values = SHP_SET, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  labs(title = "Within-T21 cis-eQTL: effect size and significance per gene, with controls",
       subtitle = paste0(
         "Slope of expression on alt-allele dosage at the best variant (sign follows the alt allele, which is arbitrary) against the gene-level\n",
         "permutation q. With 1000 permutations q cannot fall much below 0.001, so the strongest genes share the top row.",
         if (length(untested)) paste0("\nNo GTEx variants, not tested: ", paste(untested, collapse = ", "), ".") else ""),
       x = "Within-T21 slope at the best variant (log2-CPM per alt allele)",
       y = expression(-log[10]~italic(q)~"(gene-level permutation test)")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 7.3, colour = "grey30"),
        legend.position = "bottom", legend.text = element_text(size = 7.5))
save_fig(fig_effects, run$figure("eqtl_effect_sizes"), width = 8.5, height = 6.8)

# =============================================================================
# STEP 4: chr21 map
# =============================================================================

cat("\nStep 4: chr21 map...\n")
bands <- map_bands(lanes, roster, extra)
bands[, tss_mb := tss / 1e6]
bands[, outcome := factor(LAB_EQTL[eqtl_lane], levels = LAB_EQTL)]
bands[, direction := factor(direction, levels = c("up", "down"))]
print(bands[, .(Gene_name, tss_mb = round(tss_mb, 2), direction, eqtl_lane)])

chrom <- data.table(xmin = c(0, CEN[1], CEN[2]) / 1e6,
                    xmax = c(CEN[1], CEN[2], CHR21_LEN) / 1e6,
                    part = c("p arm", "centromere", "q arm"))
BAND_W <- 0.12   # Mb, half-width of a gene band
up   <- bands[direction == "up"]
down <- bands[direction == "down"]

fig_map <- ggplot() +
  geom_rect(data = chrom, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1),
            fill = c("grey85", "grey55", "grey92"), colour = "grey40", linewidth = 0.3) +
  geom_rect(data = bands, aes(xmin = tss_mb - BAND_W, xmax = tss_mb + BAND_W,
                              ymin = 0, ymax = 1, fill = outcome), colour = NA) +
  geom_segment(data = up,   aes(x = tss_mb, xend = tss_mb, y = 1, yend = 1.25), colour = "grey50", linewidth = 0.3) +
  geom_segment(data = down, aes(x = tss_mb, xend = tss_mb, y = 0, yend = -0.25), colour = "grey50", linewidth = 0.3) +
  geom_text_repel(data = up, aes(x = tss_mb, y = 1.3, label = Gene_name, colour = direction),
                  direction = "x", angle = 90, hjust = 0, vjust = 0.5, size = 2.6,
                  segment.size = 0.25, segment.colour = "grey60", box.padding = 0.15,
                  max.overlaps = Inf, seed = 1, ylim = c(1.3, NA), show.legend = FALSE) +
  geom_text_repel(data = down, aes(x = tss_mb, y = -0.3, label = Gene_name, colour = direction),
                  direction = "x", angle = 90, hjust = 1, vjust = 0.5, size = 2.6,
                  segment.size = 0.25, segment.colour = "grey60", box.padding = 0.15,
                  max.overlaps = Inf, seed = 1, ylim = c(NA, -0.3), show.legend = FALSE) +
  # legend keys for direction, drawn from a zero-size layer
  geom_point(data = data.frame(x = c(-10, -10), y = c(0.5, 0.5), direction = factor(c("up", "down"), levels = c("up", "down"))),
             aes(x = x, y = y, colour = direction), size = 0) +
  scale_fill_manual(values = setNames(COL_EQTL, LAB_EQTL), drop = FALSE, name = "Gene-level cis-eQTL test") +
  scale_colour_manual(values = COL_DIR, labels = c(up = "Higher than expected (DE_high)", down = "Lower than expected (DE_low)"),
                      name = "Deviation", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15)),
         fill = guide_legend(override.aes = list(colour = NA))) +
  scale_x_continuous(breaks = seq(0, 45, 5), expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = c(0, CHR21_LEN / 1e6), ylim = c(-3.2, 4.2), clip = "off") +
  labs(title = "Deviating chr21 genes and their cis-eQTL outcome",
       subtitle = sprintf("Bands at the TSS of each of the %d deviating genes (GRCh38); labels above the bar are higher, below are lower. Centromere in dark grey.", nrow(bands)),
       x = "chr21 position (Mb)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        axis.line.x = element_line(colour = "grey40"), axis.ticks.x = element_line(colour = "grey40"),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 7.5, colour = "grey30"),
        legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 7.5))
save_fig(fig_map, run$figure("chr21_deviating_map"), width = 10, height = 4.2)

writeLines(capture.output(sessionInfo()), run$figure("eqtl_figures_session_info.txt"))
cat("\n=== eQTL figures complete ===\n")
