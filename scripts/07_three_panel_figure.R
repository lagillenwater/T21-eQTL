# 07_three_panel_figure.R
#
# Purpose: Chr21_DEG - a 2 x 2 volcano figure.
#   A  uncorrected DESeq2, all protein-coding genes, chr21 highlighted
#   B  ploidy-corrected DESeq2, all protein-coding genes, chr21 highlighted
#   C  uncorrected, chr21 genes only
#   D  ploidy-corrected, chr21 genes only
#
# All four panels share axes so the shift is readable directly: before
# correction the chr21 cloud sits near log2(1.5) = 0.585 and is almost
# uniformly significant; after correction it recenters on 0 and most of that
# significance is absorbed, leaving the genes that genuinely deviate from the
# trisomy expectation. The bottom row removes the genome-wide background so the
# chr21 structure is visible on its own. Labelled genes are read from the lane
# table rather than hardcoded, so the figure cannot drift from the pipeline.
#
# The lane-flow diagram is NOT drawn here: paste
# results/tables/chr21_lane_sankeymatic_input.txt into
# https://sankeymatic.com/build/ to render it.
#
# Inputs:
#   - results/tables/deseq2_all_genes_both_analyses.csv (script 01)
#   - results/tables/chr21_lane_assignments.csv         (script 04)
# Outputs:
#   - results/figures/Chr21_DEG.pdf
#   - results/figures/Chr21_DEG.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---- constants --------------------------------------------------------------
ALPHA       <- 0.01         # padj threshold, matches scripts 02/04
TRISOMY_LFC <- log2(1.5)    # 0.585, the expected chr21 dosage bump
Y_CAP       <- 60           # -log10(padj) display ceiling; see note below
OUT_STEM    <- "results/figures/Chr21_DEG"

# Group labels, ordered high -> low so the legend reads top-down. Every
# (sig_lane, eqtl_lane) combination the pipeline can produce needs a label
# here: a combination with no branch falls through to NA and is silently drawn
# as "Other chr21" grey. LAB_HI_UNEX is empty at present but is the
# OLIG2-shaped case on the high side, and would otherwise disappear.
LAB_HI_EXPL <- "Higher, cis-eQTL detected"
LAB_HI_UNEX <- "Higher, no cis-eQTL detected"
LAB_HI_NOEQ <- "Higher, no GTEx eQTL"
LAB_LO_EXPL <- "Lower, cis-eQTL detected"
LAB_LO_UNEX <- "Lower, no cis-eQTL detected"
LAB_LO_NOEQ <- "Lower, no GTEx eQTL"
LAB_NOPADJ  <- "Outside noise, padj not estimable"
LAB_CHR21   <- "Other chr21"
LAB_OTHER   <- "Other protein-coding"

PALETTE <- c("#B2182B", "#67001F", "#E08214", "#2166AC", "#D6604D", "#762A83",
             "#1B7837", "#F4A582", "grey80")
names(PALETTE) <- c(LAB_HI_EXPL, LAB_HI_UNEX, LAB_HI_NOEQ, LAB_LO_EXPL,
                    LAB_LO_UNEX, LAB_LO_NOEQ, LAB_NOPADJ, LAB_CHR21, LAB_OTHER)

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
    sig_lane == "DE_high" & eqtl_lane == "cis_eqtl"     ~ LAB_HI_EXPL,
    sig_lane == "DE_high" & eqtl_lane == "no_cis_eqtl"  ~ LAB_HI_UNEX,
    sig_lane == "DE_high" & eqtl_lane == "no_GTEx_data" ~ LAB_HI_NOEQ,
    sig_lane == "DE_low"  & eqtl_lane == "cis_eqtl"     ~ LAB_LO_EXPL,
    sig_lane == "DE_low"  & eqtl_lane == "no_cis_eqtl"  ~ LAB_LO_UNEX,
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

# ---- panels C and D: chr21 only -------------------------------------------
# Same builder, same axes, but the genome-wide background is removed so the
# 160 chr21 genes can be read on their own.
res_chr21 <- res %>% filter(Chr == "chr21")

panel_c <- build_volcano(
  res_chr21, "raw_log2FC", "raw_padj",
  title    = "C  Uncorrected, chr21 only",
  subtitle = sprintf("%d protein-coding chr21 genes; dotted = 1.5x expectation",
                     nrow(res_chr21)),
  show_trisomy_line = TRUE
)

panel_d <- build_volcano(
  res_chr21, "norm_log2FC", "norm_padj",
  title    = "D  Ploidy-corrected, chr21 only",
  subtitle = sprintf("%d protein-coding chr21 genes; deviating = %s",
                     nrow(res_chr21),
                     paste(sort(lane$Gene_name[lane$sig_lane %in% c("DE_high", "DE_low")]),
                           collapse = ", ")),
  show_trisomy_line = FALSE
)

# ---- assemble ---------------------------------------------------------------
# One legend, deterministically: patchwork only merges guides it considers
# identical, and the four panels build theirs from different data. Keep the
# legend on panel A only and let guides = "collect" place it in guide_area().
panel_a <- panel_a + theme(legend.position = "bottom")
for (nm in c("panel_b", "panel_c", "panel_d")) {
  assign(nm, get(nm) + theme(legend.position = "none"))
}

design <- "AB
CD
EE"

fig <- panel_a + panel_b + panel_c + panel_d + guide_area() +
  plot_layout(design = design, guides = "collect",
              heights = c(1, 1, 0.14))

cat("Writing figure...\n")
# Plain pdf() rather than cairo_pdf: cairo is not available on every machine
# here (no X11), and cairo_pdf fails to write at all when it is missing.
ggsave(paste0(OUT_STEM, ".pdf"), fig, width = 9.5, height = 9.5,
       units = "in", device = "pdf")
ggsave(paste0(OUT_STEM, ".png"), fig, width = 9.5, height = 9.5,
       units = "in", dpi = 300)

# Confirm both actually landed - a failed graphics device is otherwise silent.
for (f in paste0(OUT_STEM, c(".pdf", ".png"))) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat(sprintf("  Saved: %s (%.1f KB)\n", f, file.size(f) / 1024))
}

writeLines(capture.output(sessionInfo()),
           "results/figures/Chr21_DEG_session_info.txt")
cat("Done.\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  ADDED LAB_HI_UNEX ("Higher, no cis-eQTL detected") with its own
#             palette entry and case_when branch. DE_high & no_cis_eqtl had no
#             branch, so such a gene fell through to NA and was drawn as grey
#             "Other chr21" - silently unlabelled. The combination is empty
#             today but is the OLIG2-shaped case on the high side, and the two
#             DE_high genes are one permutation result away from it.
#
# 2026-08-31  RELABELLED the panel C sub-categories: "DE high (>= 1.5 raw FC)"
#             / "DE low (< 1.5 raw FC)" named the raw fold change, but the rule
#             is abs(norm_log2FC) >= log2(1.5) on the ploidy-corrected scale.
#             The panel title, subtitle and header comment likewise described
#             the retired cohort-noise SD filter and a stale gene count; they
#             now state the current rule without hardcoded numbers.
# 2026-09-01  REPLACED panel C (the ggalluvial lane-flow) with two chr21-only
#             volcano panels (C uncorrected, D ploidy-corrected), making a 2x2
#             figure. Output renamed three_panel_summary -> Chr21_DEG. The lane
#             flow is rendered by hand from chr21_lane_sankeymatic_input.txt at
#             sankeymatic.com instead. Dropped the now-unused ggalluvial, png,
#             grid and tidyr loads and the SANKEY_PNG switch.
#             Reason: user request - the SankeyMATIC render is preferred and
#             the chr21-only view is more legible than the flattened alluvial.
