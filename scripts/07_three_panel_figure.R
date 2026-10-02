# 07_three_panel_figure.R
#
# Purpose: Two volcano figures from the DESeq2 results.
#   volcano_all_genes  A uncorrected / B ploidy-corrected, every gene of the
#                      target biotypes in grey with chr21 genes drawn on top
#                      in one highlight colour. No labels: the point is the
#                      shift of the whole chr21 cloud from log2(1.5) to 0.
#   volcano_chr21      A uncorrected / B ploidy-corrected, chr21 genes only.
#                      Deviating genes are labelled and coloured by direction
#                      alone (red DE_high, blue DE_low; tier 1 bold); every
#                      other chr21 gene is grey. eQTL outcomes are not drawn
#                      here (see scripts/11_eqtl_figures.R).
#
# Both figures share axes so the shift is readable directly: before correction
# the chr21 cloud sits near log2(1.5) = 0.585 and is almost uniformly
# significant; after correction it recenters on 0 and most of that
# significance is absorbed, leaving the genes that deviate from the trisomy
# expectation. Labelled genes are read from the lane table, not hardcoded.
#
# The lane-flow diagram is NOT drawn here: paste
# results/runs/<name>/tables/chr21_lane_sankeymatic_input.txt into
# https://sankeymatic.com/build/ to render it.
#
# Inputs (per run):
#   - tables/deseq2_all_genes_both_analyses.csv (script 01)
#   - tables/chr21_lane_assignments.csv         (script 04)
# Outputs (per run):
#   - figures/volcano_all_genes.{pdf,png}
#   - figures/volcano_chr21.{pdf,png}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---- constants --------------------------------------------------------------
source("scripts/lib/run.R"); run <- load_run()
ALPHA       <- run$thresholds$alpha_de   # padj threshold, from the run config
TRISOMY_LFC <- log2(run$ploidy)          # 0.585, the expected chr21 dosage bump
Y_CAP       <- 60           # -log10(padj) display ceiling; see note below
source("scripts/lib/biotypes.R")   # TARGET_BIOTYPES, matches scripts 02/04/06

LAB_HIGH  <- "Higher than expected (DE_high)"
LAB_LOW   <- "Lower than expected (DE_low)"
LAB_CHR21 <- "Other chr21"
LAB_OTHER <- "Other genes (target biotypes)"
# Direction colours match scripts/11_eqtl_figures.R.
PALETTE <- c("#B2182B", "#2166AC", "#762A83", "grey80")
names(PALETTE) <- c(LAB_HIGH, LAB_LOW, LAB_CHR21, LAB_OTHER)
CHR21_GREY <- "grey60"   # chr21 background in the chr21-only figure

# ---- load -------------------------------------------------------------------
cat("Loading DESeq2 results and lane assignments...\n")
res  <- read_csv(run$table("deseq2_all_genes_both_analyses.csv"),
                 show_col_types = FALSE)
lane <- read_csv(run$table("chr21_lane_assignments.csv"),
                 show_col_types = FALSE)

res <- res %>% filter(Gene_type %in% TARGET_BIOTYPES)

# No global padj filter here: build_volcano() filters on its own padj column
# per panel, so a gene missing raw_padj but carrying a valid norm_padj still
# appears in the corrected panels (and vice versa).
cat(sprintf("  %d target-biotype genes (%d on chr21)\n",
            nrow(res), sum(res$Chr == "chr21")))

# ---- assign display groups from the lane table ------------------------------
# Derived, not hardcoded: whatever scripts 02/04 currently call DE is what gets
# labelled here.
lane_groups <- lane %>%
  mutate(group = case_when(sig_lane == "DE_high" ~ LAB_HIGH,
                           sig_lane == "DE_low"  ~ LAB_LOW,
                           TRUE ~ NA_character_)) %>%
  filter(!is.na(group)) %>%
  select(Gene_name, group, tier)
cat(sprintf("  Deviating genes: %d higher, %d lower (tier 1: %d)\n",
            sum(lane_groups$group == LAB_HIGH), sum(lane_groups$group == LAB_LOW),
            sum(lane_groups$tier == 1L, na.rm = TRUE)))

res <- res %>%
  left_join(lane_groups, by = "Gene_name") %>%
  mutate(group = case_when(
    !is.na(group)  ~ group,
    Chr == "chr21" ~ LAB_CHR21,
    TRUE           ~ LAB_OTHER
  ),
  group = factor(group, levels = names(PALETTE)))

# ---- volcano builder --------------------------------------------------------
# A handful of genes have padj ~1e-216 (uncorrected chr21 genes tested against
# FC = 1), which would flatten every other point against the x-axis. Cap the
# display at Y_CAP and draw capped points as triangles, keyed in the shape
# legend and counted in the panel subtitle, so the truncation is visible
# rather than silent.
#
# mode = "all":   background grey, chr21 (deviating or not) in the chr21
#                 highlight colour, no labels.
# mode = "chr21": chr21 background grey, deviating genes coloured by direction
#                 and labelled (tier 1 bold).
build_volcano <- function(df, lfc_col, padj_col, title, subtitle,
                          show_trisomy_line, mode = c("all", "chr21")) {
  mode <- match.arg(mode)
  d <- df %>%
    filter(!is.na(.data[[padj_col]])) %>%
    mutate(
      lfc      = .data[[lfc_col]],
      neglog10 = -log10(.data[[padj_col]] + 1e-300),
      capped   = factor(neglog10 > Y_CAP, levels = c("FALSE", "TRUE")),
      y        = pmin(neglog10, Y_CAP)
    )
  n_capped <- sum(d$capped == "TRUE")
  if (n_capped > 0) {
    subtitle <- sprintf("%s; %d at the ceiling", subtitle, n_capped)
  }

  if (mode == "all") {
    d <- d %>% mutate(group = factor(ifelse(Chr == "chr21", LAB_CHR21, LAB_OTHER),
                                     levels = c(LAB_CHR21, LAB_OTHER)))
    pal    <- PALETTE[c(LAB_CHR21, LAB_OTHER)]
    bg     <- d %>% filter(group == LAB_OTHER)
    c21    <- d %>% filter(group == LAB_CHR21)
    named  <- d[0, ]
    c21_col <- NULL
  } else {
    d <- d %>% filter(Chr == "chr21") %>%
      mutate(group = factor(as.character(group), levels = c(LAB_HIGH, LAB_LOW, LAB_CHR21)))
    pal    <- c(PALETTE[c(LAB_HIGH, LAB_LOW)], setNames(CHR21_GREY, LAB_CHR21))
    bg     <- d[0, ]
    c21    <- d %>% filter(group == LAB_CHR21)
    named  <- d %>% filter(group != LAB_CHR21)
  }

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
    geom_text_repel(data = named,
                    aes(label = Gene_name,
                        fontface = ifelse(!is.na(tier) & tier == 1L, "bold", "plain")),
                    size = 2.4,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.2, segment.alpha = 0.6,
                    box.padding = 0.4, force = 6, seed = 1,
                    show.legend = FALSE) +
    # limits pins the legend to the full level set so the two panels build the
    # same guide and patchwork collects one legend.
    scale_colour_manual(values = pal, limits = names(pal), drop = FALSE, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), drop = FALSE,
                       name = NULL,
                       labels = c(sprintf("-log10 padj <= %d", Y_CAP),
                                  sprintf("-log10 padj > %d, drawn at the ceiling",
                                          Y_CAP))) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1),
                                 nrow = 2, byrow = TRUE, order = 1),
           shape  = guide_legend(override.aes = list(size = 2.2, alpha = 1,
                                                     colour = "grey30"),
                                 nrow = 2, order = 2)) +
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

# One legend per figure: keep it on panel A, collect with patchwork.
assemble <- function(pa, pb) {
  # Panel B carries no legend, so guides = "collect" gathers panel A's alone;
  # setting legend.position through `&` would switch B's back on and draw
  # the shape legend twice.
  pa <- pa + theme(legend.position = "bottom")
  pb <- pb + theme(legend.position = "none")
  pa + pb + guide_area() +
    plot_layout(design = "AB\nCC", guides = "collect", heights = c(1, 0.16)) &
    theme(legend.box = "horizontal", legend.box.just = "top",
          legend.spacing.x = unit(1.5, "lines"))
}

write_fig <- function(fig, stem, width, height) {
  # Plain pdf() rather than cairo_pdf: cairo is not available on every machine
  # here (no X11), and cairo_pdf fails to write at all when it is missing.
  ggsave(paste0(stem, ".pdf"), fig, width = width, height = height, units = "in", device = "pdf")
  ggsave(paste0(stem, ".png"), fig, width = width, height = height, units = "in", dpi = 300)
  # Confirm both actually landed - a failed graphics device is otherwise silent.
  for (f in paste0(stem, c(".pdf", ".png"))) {
    if (!file.exists(f)) stop("failed to write ", f)
    cat(sprintf("  Saved: %s (%.1f KB)\n", f, file.size(f) / 1024))
  }
}

# ---- figure 1: all genes, chr21 highlighted, no labels ---------------------
fig_all <- assemble(
  build_volcano(res, "raw_log2FC", "raw_padj", mode = "all",
    title    = "A  Uncorrected",
    subtitle = sprintf("chr21 median FC = %.2f (dotted = 1.5x expectation); %d/%d chr21 padj < %.2g",
                       2^med_raw, sig_raw, n_chr21, ALPHA),
    show_trisomy_line = TRUE),
  build_volcano(res, "norm_log2FC", "norm_padj", mode = "all",
    title    = "B  Ploidy-corrected",
    subtitle = sprintf("chr21 median log2FC = %.2f; %d/%d chr21 padj < %.2g",
                       med_norm, sig_norm, n_chr21, ALPHA),
    show_trisomy_line = FALSE))
cat("Writing volcano_all_genes...\n")
write_fig(fig_all, run$figure("volcano_all_genes"), width = 9.5, height = 5.8)

# ---- figure 2: chr21 only, deviating genes labelled by direction -----------
tier1 <- sort(lane$Gene_name[!is.na(lane$tier) & lane$tier == 1L])
fig_chr21 <- assemble(
  build_volcano(res, "raw_log2FC", "raw_padj", mode = "chr21",
    title    = "A  Uncorrected, chr21 only",
    subtitle = sprintf("%d chr21 genes (%s); dotted = 1.5x expectation",
                       n_chr21, TARGET_BIOTYPES_LABEL),
    show_trisomy_line = TRUE),
  build_volcano(res, "norm_log2FC", "norm_padj", mode = "chr21",
    title    = "B  Ploidy-corrected, chr21 only",
    subtitle = sprintf("%d deviating genes labelled; tier 1 (bold):\n%s",
                       nrow(lane_groups),
                       paste(strwrap(paste(tier1, collapse = ", "), width = 95), collapse = "\n")),
    show_trisomy_line = FALSE))
cat("Writing volcano_chr21...\n")
write_fig(fig_chr21, run$figure("volcano_chr21"), width = 9.5, height = 5.8)

writeLines(capture.output(sessionInfo()),
           run$figure("volcano_session_info.txt"))
cat("Done.\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-09-10  SPLIT the 2x2 Chr21_DEG figure into volcano_all_genes (chr21
#             highlighted, no labels) and volcano_chr21 (deviating genes
#             labelled, coloured by direction only). eQTL outcomes moved to
#             scripts/11_eqtl_figures.R. Chr21_DEG.{pdf,png} are no longer
#             written.
