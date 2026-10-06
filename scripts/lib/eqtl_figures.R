# eqtl_figures.R
#
# Data preparation for scripts/11_eqtl_figures.R. Pure table-to-table
# functions so the figure script only draws:
#
#   panel_variants  which variant each per-gene dosage panel shows, by group
#   map_bands       TSS, direction and eQTL outcome per deviating gene
#   orient_panels   which allele each panel's x-axis carries, and how it
#                   relates to the canonical GTEx direction
#   panel_dosage    dosage of that allele, from ALT dosage
#   best_variant_effects  within-T21 slope, SE and p at each panel's variant
#                         (on whichever dosage coding the caller passes)
#
# Inputs are the tables scripts 03 and 04 already write (eqtl_gene_level_perm,
# chr21_lane_assignments); no statistic is recomputed here.

PANEL_GROUPS <- c("DE high", "DE low")

lane_block <- function(sig_lane) {
  data.table::fcase(sig_lane == "DE_high", "DE high",
                    sig_lane == "DE_low",  "DE low",
                    default = NA_character_)
}

#' One row per dosage panel: the deviating gene whose expression is drawn,
#' its best variant (the one the gene-level test rests on) and the block the
#' panel belongs to. Genes with no tested variant get no panel.
panel_variants <- function(perm, lanes) {
  lanes <- data.table::as.data.table(lanes)[, .(Gene_name, sig_lane)]
  out <- merge(data.table::as.data.table(perm)[!is.na(best_variant)],
               lanes, by = "Gene_name")
  out <- out[, .(Gene_name, panel_group = lane_block(sig_lane),
                 variant_id = best_variant, q_gene_bh,
                 detected = cis_eqtl_detected)]
  out <- out[!is.na(panel_group)]
  out[, panel_group := factor(panel_group, levels = PANEL_GROUPS)]
  data.table::setorder(out, panel_group, q_gene_bh, Gene_name)
  out[]
}

#' TSS, direction and eQTL outcome for every deviating gene. Positions come
#' from the roster (GTEx-derived TSS) first and from `extra` (a tracked
#' lookup for genes GTEx does not carry) second; any gene still without a
#' position is an error naming the genes.
map_bands <- function(lanes, roster, extra) {
  lanes <- data.table::as.data.table(lanes)
  dev <- lanes[sig_lane %in% c("DE_high", "DE_low"),
               .(Gene_name, direction = data.table::fifelse(sig_lane == "DE_high", "up", "down"),
                 eqtl_lane)]
  r <- data.table::as.data.table(roster)[!is.na(tss), .(Gene_name, tss_roster = as.numeric(tss))]
  e <- data.table::as.data.table(extra)[, .(Gene_name, tss_extra = as.numeric(tss))]
  dev <- merge(dev, unique(r, by = "Gene_name"), by = "Gene_name", all.x = TRUE)
  dev <- merge(dev, unique(e, by = "Gene_name"), by = "Gene_name", all.x = TRUE)
  dev[, tss := data.table::fifelse(is.na(tss_roster), tss_extra, tss_roster)]
  missing <- sort(dev[is.na(tss), Gene_name])
  if (length(missing))
    stop("no TSS for deviating gene(s): ", paste(missing, collapse = ", "),
         " - add them to data/chr21_gene_positions.csv", call. = FALSE)
  out <- dev[, .(Gene_name, tss, direction, eqtl_lane)]
  data.table::setorder(out, tss)
  out[]
}

#' Which allele each panel puts on its x-axis, and how the panel relates to
#' the canonical GTEx direction.
#'
#' A deviating-gene panel is oriented to the DEVIATION-MATCHING allele: the
#' allele whose GTEx effect has the same sign as the gene's own deviation. A
#' DE high panel then trends up and a DE low panel trends down whenever the
#' within-T21 fit reproduces GTEx, so the trend direction reads against the
#' block label instead of looking arbitrary. The orientation is taken from
#' GTEx, never from the within-T21 slope - orienting on the data being
#' plotted would force every panel to trend the right way by construction and
#' show nothing.
#'
#' The minor allele is plotted exactly when it is itself the deviation-
#' matching one, so `plot_allele_role` doubles as the with/against-the-
#' deviation badge. A panel with nothing to match (no GTEx slope for the gene
#' at that variant) stays on minor-allele dosage.
#'
#' @param panels  output of panel_variants()
#' @param variants variant-level alleles: variant_id, minor_allele,
#'   major_allele, alt_is_minor, gtex_maf
#' @param gtex gene-specific GTEx slope per minor allele: variant_id,
#'   Gene_name, gtex_slope_minor
#' @param gene_dir Gene_name, deviation_sign (+1 / -1) for deviating genes
orient_panels <- function(panels, variants, gtex, gene_dir) {
  out <- data.table::as.data.table(panels)
  out <- merge(out, unique(data.table::as.data.table(variants), by = "variant_id"),
               by = "variant_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, data.table::as.data.table(gtex), by = c("variant_id", "Gene_name"),
               all.x = TRUE, sort = FALSE)
  out <- merge(out, data.table::as.data.table(gene_dir), by = "Gene_name",
               all.x = TRUE, sort = FALSE)
  out[, canonical_dir := data.table::fcase(is.na(gtex_slope_minor), NA_character_,
                                           gtex_slope_minor > 0, "raises",
                                           default = "lowers")]
  # Does the MINOR allele's GTEx effect run the same way as the gene deviates?
  out[, with_deviation := !is.na(gtex_slope_minor) & !is.na(deviation_sign) &
                          sign(gtex_slope_minor) == deviation_sign]
  # Plot the minor allele when it is the deviation-matching one, the major
  # allele when it is not, and the minor allele wherever there is nothing to
  # match (no GTEx slope for the gene at that variant).
  out[, plot_minor := is.na(deviation_sign) | is.na(gtex_slope_minor) | with_deviation]
  out[, plot_allele := data.table::fifelse(plot_minor, minor_allele, major_allele)]
  out[, plot_allele_role := data.table::fifelse(plot_minor, "minor", "major")]
  out[is.na(deviation_sign) | is.na(gtex_slope_minor), with_deviation := NA]
  data.table::setorder(out, panel_group, q_gene_bh, Gene_name)
  out[]
}

#' Dosage of the allele a panel is oriented to, from ALT dosage.
#' `plot_minor` FALSE reflects the minor-allele dosage once more, which is the
#' major-allele count.
panel_dosage <- function(alt_dosage, alt_is_minor, plot_minor, ploidy) {
  d <- minor_dosage(alt_dosage, alt_is_minor, ploidy)
  v <- recycle_common(d, plot_minor)
  ifelse(is.na(v[[2]]) | is.na(v[[1]]), NA_real_, ifelse(v[[2]], v[[1]], ploidy - v[[1]]))
}

#' Within-T21 effect size at each panel's variant: slope, SE and p of
#' expression on dosage (fit_variants, the fit script 03 uses), one row per
#' (Gene_name, panel_group, variant_id). NA when the variant is monomorphic or
#' fewer than 3 subjects remain.
#'
#' `dosage_col` names the regressor. Script 11 passes minor-allele dosage so
#' the plotted slope is per copy of the minor allele; p, and so the check
#' against the stored test p, is the same either way, since reflecting the
#' regressor only flips the slope's sign.
#' With `covariates` (a matrix with one row per subject, rownames subject
#' ids; `long` then needs subject_id) each fit is expression ~ dosage +
#' covariates, as in script 03. `mean_dosage` and `mean_expr` place the
#' fitted line through the panel's mean point for drawing.
#' @param long data.table with Gene_name, panel_group, variant_id, expr and `dosage_col`
best_variant_effects <- function(long, dosage_col = "alt_dosage", covariates = NULL) {
  stopifnot(all(c("Gene_name", "panel_group", "variant_id", dosage_col, "expr")
                %in% names(long)),
            is.null(covariates) || "subject_id" %in% names(long))
  d <- data.table::copy(data.table::as.data.table(long))
  data.table::setnames(d, dosage_col, ".dosage")
  d <- d[!is.na(.dosage) & !is.na(expr)]
  if (is.null(d$subject_id)) d[, subject_id := NA_character_]
  d[, {
    if (.N < 3 || stats::var(.dosage) == 0) {
      list(n = .N, slope = NA_real_, se = NA_real_, p = NA_real_,
           mean_dosage = mean(.dosage), mean_expr = mean(expr))
    } else {
      C <- if (is.null(covariates)) NULL else covariates[subject_id, , drop = FALSE]
      f <- fit_variants(matrix(.dosage, ncol = 1), expr, covariates = C)
      list(n = .N, slope = f$slope[1], se = f$se[1], p = f$p[1],
           mean_dosage = mean(.dosage), mean_expr = mean(expr))
    }
  }, by = .(Gene_name, panel_group = as.character(panel_group), variant_id)]
}

#' Vertical slices for map bands that would overlap.
#'
#' Bands are drawn half_width either side of a gene's position; genes closer
#' than one band width (chained) form a cluster, and a cluster splits the bar
#' height into equal slices, bottom to top in position order, so no band hides
#' another. A lone band keeps the full height.
#' @return data.table(ymin, ymax) in the order of `pos`, on a 0-1 bar.
stack_bands <- function(pos, half_width) {
  o <- order(pos)
  p <- pos[o]
  cluster <- cumsum(c(TRUE, diff(p) >= 2 * half_width))
  k <- stats::ave(seq_along(p), cluster, FUN = seq_along)
  n <- stats::ave(seq_along(p), cluster, FUN = length)
  data.table::data.table(ymin = (k - 1) / n, ymax = k / n)[order(o)]
}

# The map's eQTL outcome. Genes without a variant to test are split by
# whether GTEx tested them in whole blood at all.
MAP_OUTCOMES <- c("cis-eQTL detected", "tested, not detected",
                  "GTEx tested, no eQTL", "not tested by GTEx")
# Filled = detected, open = tested and not detected, greys = no variant to
# test (the figure palette: red and blue are reserved for direction).
MAP_FILL <- setNames(c("grey15", "white", "grey55", "grey85"), MAP_OUTCOMES)

#' eQTL outcome of each deviating gene for the map.
#' @param eqtl_lane cis_eqtl / no_cis_eqtl / no_GTEx_data (script 04)
#' @param gtex_tested whether GTEx tested the gene in whole blood (eGenes file)
map_outcome <- function(eqtl_lane, gtex_tested) {
  bad <- setdiff(unique(eqtl_lane), c("cis_eqtl", "no_cis_eqtl", "no_GTEx_data"))
  if (length(bad)) stop("unknown eqtl_lane: ", paste(bad, collapse = ", "), call. = FALSE)
  out <- data.table::fcase(eqtl_lane == "cis_eqtl",    MAP_OUTCOMES[1],
                           eqtl_lane == "no_cis_eqtl", MAP_OUTCOMES[2],
                           gtex_tested,                MAP_OUTCOMES[3],
                           default =                   MAP_OUTCOMES[4])
  factor(out, levels = MAP_OUTCOMES)
}

#' Schematic chr21 with a band at each deviating gene's TSS (scripts 11, 12).
#'
#' Band fill = the eQTL outcome (MAP_FILL); overlapping bands share the bar
#' height (stack_bands). Gene labels, italic, are coloured by direction and
#' sit above the bar for genes higher than expected, below for lower. No
#' title: the caption carries it. Callers attach ggplot2 and ggrepel.
#' @param bands Gene_name, tss (bp), direction ("up"/"down"), outcome (map_outcome)
#' @param col_dir colours for c(up, down)
#' @param chr_len, centromere chr21 length and centromere (bp, GRCh38)
#' @param band_w half-width of a band, Mb
deviating_map_plot <- function(bands, col_dir, chr_len, centromere, band_w = 0.12) {
  b <- data.table::copy(data.table::as.data.table(bands))
  b[, tss_mb := tss / 1e6]
  b[, c("ymin", "ymax") := stack_bands(tss_mb, band_w)]
  b[, direction := factor(direction, levels = c("up", "down"))]
  up <- b[direction == "up"]; down <- b[direction == "down"]
  chrom <- data.table::data.table(xmin = c(0, centromere[1], centromere[2]) / 1e6,
                                  xmax = c(centromere[1], centromere[2], chr_len) / 1e6)
  key <- data.frame(x = c(-10, -10), y = c(0.5, 0.5), direction = factor(c("up", "down"), levels = c("up", "down")))
  ggplot() +
    geom_rect(data = chrom, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1),
              fill = c("grey85", "grey55", "grey95"), colour = "grey40", linewidth = 0.3) +
    geom_rect(data = b, aes(xmin = tss_mb - band_w, xmax = tss_mb + band_w, ymin = ymin, ymax = ymax,
                            fill = outcome), colour = "grey15", linewidth = 0.25) +
    geom_segment(data = up,   aes(x = tss_mb, xend = tss_mb, y = 1, yend = 1.25), colour = "grey50", linewidth = 0.3) +
    geom_segment(data = down, aes(x = tss_mb, xend = tss_mb, y = 0, yend = -0.25), colour = "grey50", linewidth = 0.3) +
    geom_text_repel(data = up, aes(x = tss_mb, y = 1.3, label = Gene_name, colour = direction),
                    direction = "x", angle = 90, hjust = 0, vjust = 0.5, size = 2.6, fontface = "italic",
                    segment.size = 0.25, segment.colour = "grey60", box.padding = 0.15,
                    max.overlaps = Inf, seed = 1, ylim = c(1.3, NA), show.legend = FALSE) +
    geom_text_repel(data = down, aes(x = tss_mb, y = -0.3, label = Gene_name, colour = direction),
                    direction = "x", angle = 90, hjust = 1, vjust = 0.5, size = 2.6, fontface = "italic",
                    segment.size = 0.25, segment.colour = "grey60", box.padding = 0.15,
                    max.overlaps = Inf, seed = 1, ylim = c(NA, -0.3), show.legend = FALSE) +
    # legend keys for direction, from a zero-size layer outside the bar
    geom_point(data = key, aes(x = x, y = y, colour = direction), size = 0) +
    scale_fill_manual(values = MAP_FILL, drop = FALSE, name = "Gene-level cis-eQTL test") +
    scale_colour_manual(values = col_dir, labels = c(up = "Higher than expected", down = "Lower than expected"),
                        name = "Deviation", drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3, shape = 15), order = 1),
           fill = guide_legend(override.aes = list(colour = "grey15", linewidth = 0.25), order = 2)) +
    scale_x_continuous(breaks = seq(0, 45, 5), expand = expansion(mult = c(0.01, 0.02))) +
    coord_cartesian(xlim = c(0, chr_len / 1e6), ylim = c(-3.2, 4.2), clip = "off") +
    labs(x = "chr21 position (Mb)", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          axis.line.x = element_line(colour = "grey40"), axis.ticks.x = element_line(colour = "grey40"),
          legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7.5),
          legend.spacing.y = unit(2, "pt"))
}

#' The map's genes: map_bands() plus each gene's map_outcome(), ordered by
#' position (scripts 11 and 12).
#' @param lanes lane table with ensembl_stable
#' @param gtex_tested_ids stable Ensembl ids GTEx tested in whole blood (its
#'   eGenes file lists every tested gene)
deviating_map_data <- function(lanes, roster, extra, gtex_tested_ids) {
  b <- map_bands(lanes, roster, extra)
  ids <- data.table::as.data.table(lanes)[, .(Gene_name, ensembl_stable)]
  b <- merge(b, unique(ids, by = "Gene_name"), by = "Gene_name", all.x = TRUE)
  b[, outcome := map_outcome(eqtl_lane, ensembl_stable %in% gtex_tested_ids)]
  data.table::setorder(b, tss)
  b[]
}
