# 11_eqtl_figures.R
#
# Purpose: Two figures for the within-T21 cis-eQTL stage, drawn from tables
#          scripts 03 and 04 already wrote. No statistic is recomputed here.
#
#   eqtl_dosage_panels   Expression by dosage in T21, box + jitter, one panel
#                        per gene. Main figure: panel A = DE high, panel B =
#                        DE low genes, each on the best variant of the
#                        gene-level permutation test, the one the
#                        classification rests on, so detected and
#                        not-detected genes are shown alike. Each panel's
#                        x-axis carries the DEVIATION-MATCHING allele (the
#                        allele GTEx says moves expression the way the gene
#                        deviates), so a panel reproducing GTEx trends up in
#                        panel A and down in panel B. The strip names the
#                        minor allele, its MAF, the canonical GTEx direction,
#                        whether the within-T21 fit agrees, and whether the
#                        minor allele runs with or against the deviation; a
#                        trend line coloured by that agreement is drawn on
#                        every panel.
#   chr21_deviating_map  Schematic chr21 (GRCh38) with a band at every
#                        deviating gene's TSS, coloured by eQTL outcome; gene
#                        labels coloured by direction (red up, blue down).
#
#   Also writes tables/eqtl_best_variant_effects.csv: the within-T21 slope per
#   minor allele at each gene's best variant, refit here and checked against
#   the p script 03 stored. The comparison of effect sizes with the
#   Expected-dosage genes is script 18.
#
# Inputs (per run):
#   - tables/chr21_lane_assignments.csv, eqtl_gene_level_perm.csv
#   - processed/eqtl_supported_genes.csv (TSS), genotypes_filtered.csv,
#     expression_adjusted.csv, analysis_cohort.csv
#   - data/chr21_gene_positions.csv (TSS for deviating genes GTEx lacks)
# Outputs (per run): figures/eqtl_dosage_panels.{pdf,png},
#   figures/chr21_deviating_map.{pdf,png}, tables/eqtl_best_variant_effects.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("scripts/lib/cohort.R")        # subject_id_from_labid
source("scripts/lib/eqtl_fit.R")      # expr_from_matrix
source("scripts/lib/alleles.R")       # minor_dosage, T21_CHR21_PLOIDY
source("scripts/lib/eqtl_figures.R")  # panel_variants, orient_panels, map_bands
source("scripts/lib/covariates.R")    # run_covariate_matrix: same fit as script 03
source("scripts/lib/run.R"); run <- load_run()

# Direction colours match the chr21 volcano (script 07); outcome colours are
# shared by the overview and the map. Both sets pass the colour-vision
# separation checks in light mode.
COL_DIR   <- c(up = "#B2182B", down = "#2166AC")
COL_EQTL  <- c(cis_eqtl = "#1B7837", no_cis_eqtl = "#E08214", no_GTEx_data = "grey55")
LAB_EQTL  <- c(cis_eqtl = "cis-eQTL detected", no_cis_eqtl = "Tested, not detected",
               no_GTEx_data = "No GTEx coverage")
COL_GROUP <- c("DE high" = COL_DIR[["up"]], "DE low" = COL_DIR[["down"]])

# Does the within-T21 trend reproduce the canonical (GTEx) direction of the
# eQTL for the same allele? Compared on the minor allele, so the answer is
# independent of which allele a panel is oriented to.
LAB_AGREE <- c(match    = "Within-T21 trend matches the GTEx direction",
               opposite = "Within-T21 trend opposite to GTEx")
COL_AGREE <- setNames(c("#1B7837", "#E08214"), LAB_AGREE)

# Panel-block titles and x-axis labels. A block is drawn on the allele that
# matches its deviation, so the block title states which way the panels run.
BLOCK_TITLE <- c(
  "DE high" = "DE high - expressed higher than the trisomy expectation",
  "DE low"  = "DE low - expressed lower than the trisomy expectation")
BLOCK_XLAB <- c(
  "DE high" = paste("Dosage of the expression-raising allele in T21 (0-3)",
                    "- a panel that reproduces the GTEx direction trends up"),
  "DE low"  = paste("Dosage of the expression-lowering allele in T21 (0-3)",
                    "- a panel that reproduces the GTEx direction trends down"))

# GRCh38 chr21: length and centromere (UCSC cytoBand acen bands).
CHR21_LEN <- 46709983
CEN       <- c(10864561, 12915808)

cat(sprintf("=== T21-eQTL: eQTL figures [run: %s] ===\n\n", run$name))

# =============================================================================
# STEP 1: Load tables
# =============================================================================

lanes  <- fread(run$table("chr21_lane_assignments.csv"))
perm   <- fread(run$table("eqtl_gene_level_perm.csv"))
roster <- fread(run$processed("eqtl_supported_genes.csv"))
extra  <- fread("data/chr21_gene_positions.csv")
cat(sprintf("  deviating genes: %d; tested: %d\n",
            sum(lanes$sig_lane %in% c("DE_high", "DE_low")), sum(!is.na(perm$p_gene_perm))))

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
panels <- panel_variants(perm, lanes)
cat(sprintf("  panels: %s\n",
            paste(sprintf("%s %d", levels(panels$panel_group),
                          tabulate(panels$panel_group, nbins = nlevels(panels$panel_group))),
                  collapse = "; ")))

geno     <- fread(run$processed("genotypes_filtered.csv"))
cohort   <- run$cohort()
E        <- run$expression()
meta_t21 <- cohort[Karyotype == "T21"]
meta_t21[, subject_id := subject_id_from_labid(LabID)]
# Script 03 fits expression ~ dosage + the run's covariates (projected out of
# both sides); the slopes and trend lines here use the same fit.
COV <- run_covariate_matrix(run, meta_t21)
if (!is.null(COV)) rownames(COV) <- meta_t21$subject_id

# ---- Which allele each panel is oriented to --------------------------------
# A deviating gene's panel carries the DEVIATION-MATCHING allele: the allele
# whose GTEx effect has the same sign as the gene's own deviation. Every DE
# high panel therefore trends up and every DE low panel trends down whenever
# the within-T21 fit reproduces GTEx, and a panel that does not reproduce it
# is the one trending the wrong way. The orientation comes from GTEx, never
# from the within-T21 slope being plotted.
#
# The minor allele is plotted exactly when it is itself the deviation-matching
# one, so the "minor"/"major" role printed in each strip is also the
# with/against-the-deviation badge.
tv <- fread(run$processed("eqtl_target_variants.csv"))
tv[, gtex_slope_minor := align_slope_to_minor(slope, alt_is_minor)]
variant_alleles <- unique(tv[, .(variant_id, minor_allele, major_allele,
                                 alt_is_minor, gtex_maf)], by = "variant_id")
gtex_by_gene    <- unique(tv[, .(variant_id, Gene_name, gtex_slope_minor)],
                          by = c("variant_id", "Gene_name"))
gene_dir <- lanes[sig_lane %in% c("DE_high", "DE_low"),
                  .(Gene_name, deviation_sign = sign(norm_log2FC))]

missing_allele <- setdiff(panels$variant_id, variant_alleles$variant_id)
if (length(missing_allele))
  stop("no minor-allele call for panel variant(s): ",
       paste(head(missing_allele, 3), collapse = ", "), call. = FALSE)

panels <- orient_panels(panels, variant_alleles, gtex_by_gene, gene_dir)
cat(sprintf("  minor allele runs with the gene's deviation: %d of %d deviating-gene panels\n",
            sum(panels$with_deviation, na.rm = TRUE),
            sum(!is.na(panels$with_deviation))))
cat(sprintf("  panels drawn on major-allele dosage (minor allele runs against): %d\n",
            sum(!panels$plot_minor)))

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
# dosage_minor: the issue #3 coding, what every reported slope is per.
# dosage_plot:  what this panel's x-axis shows (the same thing reflected
#               again where the major allele is the deviation-matching one).
plot_df[, dosage_minor := minor_dosage(alt_dosage, alt_is_minor, T21_CHR21_PLOIDY)]
plot_df[, dosage_plot  := panel_dosage(alt_dosage, alt_is_minor, plot_minor, T21_CHR21_PLOIDY)]

# ---- Does the within-T21 trend reproduce the canonical GTEx direction? -----
# Compared on the minor allele, so the answer does not depend on which allele
# the panel happens to be oriented to.
effects <- best_variant_effects(plot_df, dosage_col = "dosage_minor", covariates = COV)
panels[, pg := as.character(panel_group)]
panels <- merge(panels,
                effects[, .(Gene_name, pg = panel_group, variant_id, t21_slope_minor = slope)],
                by = c("Gene_name", "pg", "variant_id"), all.x = TRUE, sort = FALSE)
panels[, t21_dir := fcase(is.na(t21_slope_minor), NA_character_,
                          t21_slope_minor > 0, "raises", default = "lowers")]
panels[, agrees := canonical_dir == t21_dir]
panels[, agree_lab := fcase(is.na(agrees), NA_character_,
                            agrees,        LAB_AGREE[["match"]],
                            default =      LAB_AGREE[["opposite"]])]
cat(sprintf("  within-T21 trend reproduces the GTEx direction: %d of %d panels with one\n",
            sum(panels$agrees, na.rm = TRUE), sum(!is.na(panels$agrees))))

# ---- Facet strips -----------------------------------------------------------
setorder(panels, panel_group, q_gene_bh, Gene_name)
# Four short lines rather than two long ones: at PANEL_W inches per panel a
# strip clips past roughly 50 characters.
x_txt <- fifelse(panels$plot_allele_role == "minor",
                 sprintf("x = %s   (minor allele, MAF %.2f)",
                         panels$plot_allele, panels$gtex_maf),
                 sprintf("x = %s   (major allele; minor %s, MAF %.2f)",
                         panels$plot_allele, panels$minor_allele, panels$gtex_maf))
canon <- fifelse(is.na(panels$canonical_dir), "no GTEx direction",
                 sprintf("GTEx: %s %s", panels$minor_allele, panels$canonical_dir))
agree_txt <- fcase(is.na(panels$agrees), "", panels$agrees, "T21 agrees",
                   default = "T21 opposite")
badge <- fcase(is.na(panels$with_deviation), "",
               panels$with_deviation, "with deviation",
               default = "against deviation")
line4 <- canon
line4 <- fifelse(nzchar(agree_txt), paste0(line4, " | ", agree_txt), line4)
line4 <- fifelse(nzchar(badge),     paste0(line4, " | ", badge),     line4)
panels[, facet := sprintf("%s   q = %.2g\n%s\n%s\n%s", Gene_name, q_gene_bh,
                          sub("_b38$", "", variant_id), x_txt, line4)]

plot_df <- merge(plot_df, panels[, .(Gene_name, panel_group, variant_id, facet, agree_lab)],
                 by = c("Gene_name", "panel_group", "variant_id"))
plot_df[, facet := factor(facet, levels = panels$facet)]
plot_df[, dosage := factor(dosage_plot, levels = 0:3)]
plot_df[, agree_lab := factor(agree_lab, levels = LAB_AGREE)]

# Trend of the within-T21 fit on the allele each panel shows: script 03's
# slope (covariates projected out of both sides), drawn through the panel's
# mean point over its observed dosage range. Without covariates this is the
# OLS line of the plotted points. The factor x-axis puts dosage d at d + 1.
trend <- best_variant_effects(plot_df, dosage_col = "dosage_plot", covariates = COV)[!is.na(slope)]
trend <- merge(trend, plot_df[, .(lo = min(dosage_plot), hi = max(dosage_plot)),
                              by = .(Gene_name, panel_group = as.character(panel_group), variant_id)],
               by = c("Gene_name", "panel_group", "variant_id"))
trend <- merge(trend, unique(plot_df[, .(Gene_name, panel_group = as.character(panel_group), variant_id,
                                         facet, agree_lab)]),
               by = c("Gene_name", "panel_group", "variant_id"))
trend[, `:=`(x = lo + 1, xend = hi + 1,
             y = mean_expr + slope * (lo - mean_dosage), yend = mean_expr + slope * (hi - mean_dosage))]

# One row per category: a block's panels are laid out in a single row, so the
# figure grows sideways with the block rather than wrapping. PANEL_W is the
# width each panel needs for its four-line strip to render without clipping.
PANEL_W <- 3.1
#' One panel block. `show_legend` is TRUE for a single block per figure: with
#' a guide on every block patchwork collects one legend per block instead of
#' merging them, and the row of duplicates runs off both edges.
group_plot <- function(df, title, colour, xlab, show_legend) {
  ng <- uniqueN(df$facet)
  ggplot(df, aes(x = dosage, y = expr)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92", colour = "grey40", linewidth = 0.35) +
    geom_jitter(width = 0.18, size = 0.55, alpha = 0.45, colour = colour) +
    # Trend of the within-T21 fit on the allele this panel shows (`trend`).
    geom_segment(data = trend[panel_group == df$panel_group[1]],
                 aes(x = x, xend = xend, y = y, yend = yend, colour = agree_lab),
                 inherit.aes = FALSE, linewidth = 0.55, show.legend = show_legend) +
    scale_x_discrete(drop = FALSE) +
    scale_colour_manual(values = COL_AGREE, drop = FALSE, name = NULL) +
    # override.aes so the two outcomes absent from a given figure still get a
    # visible coloured key, which is what tells the reader what would be flagged
    guides(colour = guide_legend(override.aes = list(linewidth = 0.9))) +
    facet_wrap(~ facet, scales = "free_y", nrow = 1) +
    labs(title = title, x = xlab, y = "Expression (log2-CPM, run artifact)") +
    theme_bw(base_size = 8.5) +
    theme(strip.text = element_text(size = 5.8, lineheight = 1.05),
          plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.text = element_text(size = 7.5))
}
#' Stack the named panel_group(s) into one lettered (A, B, ...) figure and
#' save it under `stem`. Groups with no panels are silently dropped, so a
#' figure still renders if e.g. no gene lacks GTEx coverage in a given run.
dosage_fig <- function(groups, stem) {
  present <- vapply(groups, function(g) nrow(plot_df[panel_group == g]) > 0, logical(1))
  first   <- if (any(present)) groups[present][1] else NA_character_
  gp <- lapply(groups, function(g) {
    d <- plot_df[panel_group == g]
    if (nrow(d) == 0) return(NULL)
    group_plot(d, BLOCK_TITLE[[g]], COL_GROUP[[g]], BLOCK_XLAB[[g]],
               show_legend = identical(g, first))
  })
  keep <- !vapply(gp, is.null, logical(1))
  if (!any(keep)) { cat(sprintf("  Skipped %s: no panels\n", stem)); return(invisible()) }
  widest <- max(vapply(groups[keep], function(g)
    uniqueN(plot_df[panel_group == g, facet]), numeric(1)))
  fig <- wrap_plots(gp[keep], ncol = 1) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  cat(sprintf("  %s: %d block(s), widest %d panels in one row\n", stem, sum(keep), widest))
  save_fig(fig, run$figure(stem), width = 0.9 + PANEL_W * widest,
           height = max(6, 3.6 * sum(keep)))
}
dosage_fig(c("DE high", "DE low"), "eqtl_dosage_panels")

# =============================================================================
# STEP 3: Per-gene effect sizes
# =============================================================================
# Within-T21 slope of expression on MINOR-allele dosage at each gene's best
# variant, fitted in Step 2 with fit_variants(), the fit script 03 used. The
# refit p must reproduce the stored min_p_obs, which ties every effect size
# to the test it came from. Reflecting the regressor flips the slope's sign
# and nothing else, so that check holds on either coding while the sign
# means "per copy of the rarer allele" for every gene.

cat("\nStep 3: Effect sizes...\n")
chk <- merge(effects, perm[!is.na(best_variant), .(Gene_name, variant_id = best_variant, min_p_obs)],
             by = c("Gene_name", "variant_id"))
stopifnot(nrow(chk) == nrow(effects),
          max(abs(log10(chk$p) - log10(chk$min_p_obs))) < 1e-6)
cat(sprintf("  Refit p matches the stored best-variant p for all %d panels\n", nrow(chk)))
fwrite(effects, run$table("eqtl_best_variant_effects.csv"))
cat("  Wrote", run$table("eqtl_best_variant_effects.csv"), "\n")

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
