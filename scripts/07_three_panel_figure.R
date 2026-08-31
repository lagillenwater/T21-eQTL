# 07_three_panel_figure.R
#
# Purpose: Assemble the 3-panel summary figure.
#   Panel A - volcano of the UNCORRECTED (no ploidy normalization) DESeq2 results.
#   Panel B - volcano of the PLOIDY-CORRECTED DESeq2 results.
#   Panel C - flow of the chr21 genes that fall outside cohort noise, from
#             sub-category to eQTL outcome.
#
# All panels are restricted to protein-coding genes (chr21 and the genome-wide
# background alike), matching RESTRICT_TO_PROTEIN_CODING in scripts 02/04.
#
# Panels A and B share axes so the shift is readable directly: before
# correction the chr21 cloud sits near log2(1.5) = 0.585 and is almost
# uniformly significant; after correction it recenters on 0 and most of that
# significance is absorbed, leaving the genes that genuinely deviate from the
# trisomy expectation. Labelled genes are read from the lane table rather than
# hardcoded, so the figure cannot drift from the pipeline.
#
# Panel C shows only the genes OUTSIDE cohort noise. The 127 Expected-dosage
# genes are stated in the subtitle instead of drawn: at ~79% of the total they
# flatten every other lane into an illegible sliver (this is what the
# script 05 alluvial does). To use a hand-rendered SankeyMATIC diagram of the
# full flow instead, render results/tables/chr21_lane_sankeymatic_input.txt at
# https://sankeymatic.com/build/, save the PNG, and set SANKEY_PNG below.
#
# Inputs:
#   - results/tables/deseq2_all_genes_both_analyses.csv (script 01)
#   - results/tables/chr21_lane_assignments.csv         (script 04)
# Outputs:
#   - results/figures/three_panel_summary.pdf
#   - results/figures/three_panel_summary.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggalluvial)
  library(patchwork)
  library(png)
  library(grid)
})

# ---- constants --------------------------------------------------------------
ALPHA       <- 0.01         # padj threshold, matches scripts 02/04
TRISOMY_LFC <- log2(1.5)    # 0.585, the expected chr21 dosage bump
Y_CAP       <- 60           # -log10(padj) display ceiling; see note below
SANKEY_PNG  <- NULL         # set to a SankeyMATIC PNG path to use it for panel C
OUT_STEM    <- "results/figures/three_panel_summary"

# Group labels, ordered high -> low so the legend reads top-down.
LAB_HI_EXPL <- "Higher, eQTL-explained"
LAB_HI_NOEQ <- "Higher, no GTEx eQTL"
LAB_LO_EXPL <- "Lower, eQTL-explained"
LAB_LO_UNEX <- "Lower, no genetic explanation"
LAB_LO_NOEQ <- "Lower, no GTEx eQTL"
LAB_NOPADJ  <- "Outside noise, padj not estimable"
LAB_CHR21   <- "Other chr21"
LAB_OTHER   <- "Other protein-coding"

PALETTE <- c("#B2182B", "#E08214", "#2166AC", "#D6604D", "#762A83",
             "#1B7837", "#F4A582", "grey80")
names(PALETTE) <- c(LAB_HI_EXPL, LAB_HI_NOEQ, LAB_LO_EXPL, LAB_LO_UNEX,
                    LAB_LO_NOEQ, LAB_NOPADJ, LAB_CHR21, LAB_OTHER)

# ---- load -------------------------------------------------------------------
cat("Loading DESeq2 results and lane assignments...\n")
res  <- read_csv("results/tables/deseq2_all_genes_both_analyses.csv",
                 show_col_types = FALSE)
lane <- read_csv("results/tables/chr21_lane_assignments.csv",
                 show_col_types = FALSE)

res <- res %>% filter(Gene_type == "protein_coding")

# Volcano panels need a padj to place a point on the y-axis. Keep genes with a
# raw padj even when the corrected padj is NA (MX1) - dropping them silently
# would hide exactly the case worth showing.
res <- res %>% filter(!is.na(raw_padj))
cat(sprintf("  %d protein-coding genes (%d on chr21)\n",
            nrow(res), sum(res$Chr == "chr21")))

# ---- assign display groups from the lane table ------------------------------
# Derived, not hardcoded: whatever scripts 02/04 currently call DE is what gets
# labelled here.
lane_groups <- lane %>%
  mutate(group = case_when(
    sig_lane == "DE_high" & eqtl_lane == "explained"    ~ LAB_HI_EXPL,
    sig_lane == "DE_high" & eqtl_lane == "no_GTEx_data" ~ LAB_HI_NOEQ,
    sig_lane == "DE_low"  & eqtl_lane == "explained"    ~ LAB_LO_EXPL,
    sig_lane == "DE_low"  & eqtl_lane == "unexplained"  ~ LAB_LO_UNEX,
    sig_lane == "DE_low"  & eqtl_lane == "no_GTEx_data" ~ LAB_LO_NOEQ,
    sig_lane == "Not_DE_outside_noise" & is.na(norm_padj) ~ LAB_NOPADJ,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  select(Gene_name, group)

cat("  Labelled genes by group:\n")
for (g in names(PALETTE)) {
  n <- sum(lane_groups$group == g)
  if (n > 0) cat(sprintf("    %-34s %d\n", g, n))
}

res <- res %>%
  left_join(lane_groups, by = "Gene_name") %>%
  mutate(group = case_when(
    !is.na(group)  ~ group,
    Chr == "chr21" ~ LAB_CHR21,
    TRUE           ~ LAB_OTHER
  ),
  group = factor(group, levels = names(PALETTE)))

# ---- volcano builder --------------------------------------------------------
# A handful of genes have padj ~1e-216, which would flatten every other point
# against the x-axis. Cap the display at Y_CAP and draw capped points as
# triangles so the truncation is visible rather than silent.
build_volcano <- function(df, lfc_col, padj_col, title, subtitle,
                          show_trisomy_line) {
  d <- df %>%
    filter(!is.na(.data[[padj_col]])) %>%
    mutate(
      lfc      = .data[[lfc_col]],
      neglog10 = -log10(.data[[padj_col]] + 1e-300),
      capped   = neglog10 > Y_CAP,
      y        = pmin(neglog10, Y_CAP)
    )

  # Draw in layers so labelled genes sit above the chr21 cloud, which sits
  # above the genome-wide background.
  bg    <- d %>% filter(group == LAB_OTHER)
  c21   <- d %>% filter(group == LAB_CHR21)
  named <- d %>% filter(!group %in% c(LAB_OTHER, LAB_CHR21))

  ggplot(mapping = aes(x = lfc, y = y, colour = group)) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
               colour = "grey45", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey45", linewidth = 0.3) +
    {if (show_trisomy_line)
      geom_vline(xintercept = TRISOMY_LFC, linetype = "dotted",
                 colour = "#2166AC", linewidth = 0.45)} +
    geom_point(data = bg,    aes(shape = capped), size = 0.5, alpha = 0.45) +
    geom_point(data = c21,   aes(shape = capped), size = 0.9, alpha = 0.75) +
    geom_point(data = named, aes(shape = capped), size = 1.8) +
    geom_text_repel(data = named, aes(label = Gene_name), size = 2.4,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.2, segment.alpha = 0.6,
                    box.padding = 0.4, force = 6, seed = 1,
                    show.legend = FALSE) +
    # limits pins the legend to the full level set. Without it the two panels
    # build different guides (MX1 has no corrected padj, so panel B never sees
    # that level) and patchwork's guide collection draws the legend twice.
    scale_colour_manual(values = PALETTE, limits = names(PALETTE),
                        drop = FALSE, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1),
                                 nrow = 3, byrow = TRUE)) +
    # Headroom above the cap so labels on capped points have somewhere to go.
    coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, Y_CAP + 8)) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~fold~change~(T21~vs~Control)),
         y = expression(-log[10]~adjusted~italic(p))) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 10),
      plot.subtitle    = element_text(size = 7.5, colour = "grey30"),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.75, "lines")
    )
}

c21_res  <- res %>% filter(Chr == "chr21")
med_raw  <- median(c21_res$raw_log2FC, na.rm = TRUE)
med_norm <- median(c21_res$norm_log2FC, na.rm = TRUE)
sig_raw  <- sum(c21_res$raw_padj  < ALPHA, na.rm = TRUE)
sig_norm <- sum(c21_res$norm_padj < ALPHA, na.rm = TRUE)
n_chr21  <- nrow(c21_res)

cat(sprintf("  chr21 median log2FC: raw %.3f (FC %.2f) -> corrected %.3f\n",
            med_raw, 2^med_raw, med_norm))
cat(sprintf("  chr21 padj < %.2g: raw %d -> corrected %d (of %d)\n",
            ALPHA, sig_raw, sig_norm, n_chr21))

panel_a <- build_volcano(
  res, "raw_log2FC", "raw_padj",
  title    = "A  Uncorrected",
  subtitle = sprintf(
    "chr21 median FC = %.2f (dotted = 1.5x expectation); %d/%d chr21 padj < %.2g",
    2^med_raw, sig_raw, n_chr21, ALPHA),
  show_trisomy_line = TRUE
)

panel_b <- build_volcano(
  res, "norm_log2FC", "norm_padj",
  title    = "B  Ploidy-corrected",
  subtitle = sprintf("chr21 median log2FC = %.2f; %d/%d chr21 padj < %.2g",
                     med_norm, sig_norm, n_chr21, ALPHA),
  show_trisomy_line = FALSE
)

# ---- panel C ----------------------------------------------------------------
n_total    <- nrow(lane)
n_expected <- sum(lane$sig_lane == "Expected_dosage")
n_outside  <- n_total - n_expected

if (!is.null(SANKEY_PNG) && file.exists(SANKEY_PNG)) {
  cat(sprintf("Panel C: using SankeyMATIC raster %s\n", SANKEY_PNG))
  panel_c <- wrap_elements(full = grobTree(
    rasterGrob(readPNG(SANKEY_PNG), interpolate = TRUE)
  )) +
    labs(title = "C  Chr21 gene classification and eQTL support") +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0,
                                    margin = margin(b = 4)))
} else {
  cat("Panel C: drawing alluvial of the outside-cohort-noise genes\n")

  pretty_eqtl <- c(explained    = "eQTL-explained",
                   unexplained  = "eQTL-tested,\nnot supported",
                   no_GTEx_data = "No GTEx\neQTL data",
                   not_evaluated = "Not evaluated")

  flow <- lane %>%
    filter(sig_lane != "Expected_dosage") %>%
    mutate(
      Subcategory = recode(sig_lane,
                           DE_high              = "DE high\n(>= 1.5 raw FC)",
                           DE_low               = "DE low\n(< 1.5 raw FC)",
                           High_repeats         = "High repeats",
                           Low_expression       = "Low expression",
                           Not_DE_outside_noise = "Not DE"),
      Outcome = unname(pretty_eqtl[eqtl_lane])
    ) %>%
    count(Subcategory, Outcome, name = "n")

  # Order strata by size so the DE lanes, the point of the panel, sit on top.
  sub_levels <- flow %>% group_by(Subcategory) %>% summarise(t = sum(n)) %>%
    arrange(desc(t)) %>% pull(Subcategory)
  out_levels <- flow %>% group_by(Outcome) %>% summarise(t = sum(n)) %>%
    arrange(desc(t)) %>% pull(Outcome)

  flow <- flow %>%
    mutate(Subcategory = factor(Subcategory, levels = rev(sub_levels)),
           Outcome     = factor(Outcome,     levels = rev(out_levels)))

  panel_c <- ggplot(flow,
                    aes(axis1 = Subcategory, axis2 = Outcome, y = n)) +
    geom_alluvium(aes(fill = Subcategory), width = 0.28, alpha = 0.75,
                  colour = NA) +
    geom_stratum(width = 0.28, fill = "grey96", colour = "grey40",
                 linewidth = 0.3) +
    # Strata of 1-3 genes are too thin to hold text without colliding with
    # their neighbours, so label those outside the stratum with a leader line
    # and keep in-stratum text for the rest.
    geom_text(stat = "stratum", min.y = 4,
              aes(label = sprintf("%s\n(%d)", after_stat(stratum),
                                  after_stat(count))),
              size = 2.3, lineheight = 0.9) +
    geom_text_repel(stat = "stratum", max.y = 4,
                    aes(label = sprintf("%s (%d)", after_stat(stratum),
                                        after_stat(count))),
                    size = 2.2, direction = "y", nudge_x = 0.35,
                    segment.size = 0.2, segment.colour = "grey50",
                    min.segment.length = 0, box.padding = 0.15, seed = 1) +
    scale_x_discrete(limits = c("Sub-category", "eQTL outcome"),
                     expand = c(0.16, 0.16)) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = "C  Fate of the chr21 genes outside cohort noise",
      subtitle = sprintf(
        paste0("%d of %d chr21 protein-coding genes deviate by >= 1 cohort-noise SD. ",
               "The other %d are within noise (expected dosage) and not shown."),
        n_outside, n_total, n_expected),
      y = "Number of genes", x = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", size = 10),
      plot.subtitle      = element_text(size = 7.5, colour = "grey30"),
      axis.text.x        = element_text(face = "bold", size = 8)
    )
}

# ---- assemble ---------------------------------------------------------------
# One legend, deterministically. patchwork's guides="collect" only merges
# guides it considers identical, and these two are not: panel B drops MX1
# (no corrected padj), so its colour guide is built from different data and
# both legends get drawn. Rather than rely on the merge, strip panel B's
# legend outright and let the single remaining one land in guide_area().
panel_a <- panel_a + theme(legend.position = "bottom")
panel_b <- panel_b + theme(legend.position = "none")

design <- "AB
CC
DD"

fig <- panel_a + panel_b + guide_area() + panel_c +
  plot_layout(design = design, guides = "collect",
              heights = c(1, 0.14, 0.88))

cat("Writing figure...\n")
# Plain pdf() rather than cairo_pdf: cairo is not available on every machine
# here (no X11), and cairo_pdf fails to write at all when it is missing.
ggsave(paste0(OUT_STEM, ".pdf"), fig, width = 9.5, height = 10,
       units = "in", device = "pdf")
ggsave(paste0(OUT_STEM, ".png"), fig, width = 9.5, height = 10,
       units = "in", dpi = 300)

# Confirm both actually landed - a failed graphics device is otherwise silent.
for (f in paste0(OUT_STEM, c(".pdf", ".png"))) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat(sprintf("  Saved: %s (%.1f KB)\n", f, file.size(f) / 1024))
}

writeLines(capture.output(sessionInfo()),
           "results/figures/three_panel_summary_session_info.txt")
cat("Done.\n")
