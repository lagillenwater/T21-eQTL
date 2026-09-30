# eqtl_figures.R
#
# Data preparation for scripts/11_eqtl_figures.R. Pure table-to-table
# functions so the figure script only draws:
#
#   panel_variants  which variant each per-gene dosage panel shows, by group
#   overview_rows   one row per gene: own-variant q, decoy q, outcome, block
#   map_bands       TSS, direction and eQTL outcome per deviating gene
#   orient_panels   which allele each panel's x-axis carries, and how it
#                   relates to the canonical GTEx direction
#   panel_dosage    dosage of that allele, from ALT dosage
#   best_variant_effects  within-T21 slope, SE and p at each panel's variant
#                         (on whichever dosage coding the caller passes)
#   attach_effects  own-variant and decoy effects onto the overview rows
#
# Inputs are the tables scripts 03 and 04 already write (eqtl_gene_level_perm,
# eqtl_control_negative, eqtl_control_positive, chr21_lane_assignments); no
# statistic is recomputed here.

PANEL_GROUPS <- c("DE high", "DE low", "Positive control (GTEx eGenes)",
                  "Negative control (decoy variants)")
OVERVIEW_BLOCKS <- c("DE high", "DE low", "Positive control")

lane_block <- function(sig_lane) {
  data.table::fcase(sig_lane == "DE_high", "DE high",
                    sig_lane == "DE_low",  "DE low",
                    default = NA_character_)
}

#' One row per dosage panel: the gene whose expression is drawn, the variant
#' whose dosage is on the x-axis, and the group the panel belongs to.
#' Deviating genes and positive controls show their own best variant; a
#' negative-control panel shows the gene against the best variant of its
#' decoy set. Genes with no tested variant get no panel.
panel_variants <- function(perm, lanes, neg, pos) {
  lanes <- data.table::as.data.table(lanes)[, .(Gene_name, sig_lane)]
  own <- merge(data.table::as.data.table(perm)[!is.na(best_variant)],
               lanes, by = "Gene_name")
  own <- own[, .(Gene_name, panel_group = lane_block(sig_lane),
                 variant_id = best_variant, q_gene_bh,
                 detected = cis_eqtl_detected, decoy_gene = NA_character_)]
  own <- own[!is.na(panel_group)]
  p <- data.table::as.data.table(pos)[!is.na(best_variant),
         .(Gene_name, panel_group = PANEL_GROUPS[3], variant_id = best_variant,
           q_gene_bh, detected, decoy_gene = NA_character_)]
  n <- data.table::as.data.table(neg)[!is.na(best_variant),
         .(Gene_name, panel_group = PANEL_GROUPS[4], variant_id = best_variant,
           q_gene_bh, detected, decoy_gene)]
  out <- data.table::rbindlist(list(own, p, n))
  out[, panel_group := factor(panel_group, levels = PANEL_GROUPS)]
  data.table::setorder(out, panel_group, q_gene_bh, Gene_name)
  out[]
}

#' One row per gene for the overview plot: q on its own cis variants, q on
#' its decoy set, the number of variants, and the lane outcome. Deviating
#' genes come from the lane table (so genes with no GTEx coverage keep a row
#' with NA q); positive controls form their own block with no decoy.
overview_rows <- function(perm, lanes, neg, pos) {
  lanes <- data.table::as.data.table(lanes)
  dev <- lanes[sig_lane %in% c("DE_high", "DE_low"),
               .(Gene_name, block = lane_block(sig_lane), eqtl_lane)]
  dev <- merge(dev, data.table::as.data.table(perm)[, .(Gene_name, q_own = q_gene_bh,
                                                         detected_own = cis_eqtl_detected,
                                                         n_variants_own = n_variants)],
               by = "Gene_name", all.x = TRUE)
  dev <- merge(dev, data.table::as.data.table(neg)[, .(Gene_name, q_decoy = q_gene_bh,
                                                        decoy_gene)],
               by = "Gene_name", all.x = TRUE)
  dev[is.na(detected_own), detected_own := FALSE]
  p <- data.table::as.data.table(pos)[, .(Gene_name, block = OVERVIEW_BLOCKS[3],
                                          eqtl_lane = "control",
                                          q_own = q_gene_bh, detected_own = detected,
                                          n_variants_own = n_variants,
                                          q_decoy = NA_real_, decoy_gene = NA_character_)]
  out <- data.table::rbindlist(list(dev, p), use.names = TRUE)
  out[, block := factor(block, levels = OVERVIEW_BLOCKS)]
  data.table::setorder(out, block, q_own, Gene_name, na.last = TRUE)
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
#' deviation badge. Control panels (positive controls, decoy variants) have no
#' deviation to match and stay on minor-allele dosage.
#'
#' @param panels  output of panel_variants()
#' @param variants variant-level alleles: variant_id, minor_allele,
#'   major_allele, alt_is_minor, gtex_maf
#' @param gtex gene-specific GTEx slope per minor allele: variant_id,
#'   Gene_name, gtex_slope_minor. A decoy panel has no row here, because the
#'   variant is an eQTL of the decoy gene and not of the gene plotted.
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
  # match (controls, or a gene with no GTEx slope).
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
#' @param long data.table with Gene_name, panel_group, variant_id, expr and `dosage_col`
best_variant_effects <- function(long, dosage_col = "alt_dosage") {
  stopifnot(all(c("Gene_name", "panel_group", "variant_id", dosage_col, "expr")
                %in% names(long)))
  d <- data.table::copy(data.table::as.data.table(long))
  data.table::setnames(d, dosage_col, ".dosage")
  d <- d[!is.na(.dosage) & !is.na(expr)]
  d[, {
    if (.N < 3 || stats::var(.dosage) == 0) {
      list(n = .N, slope = NA_real_, se = NA_real_, p = NA_real_)
    } else {
      f <- fit_variants(matrix(.dosage, ncol = 1), expr)
      list(n = .N, slope = f$slope[1], se = f$se[1], p = f$p[1])
    }
  }, by = .(Gene_name, panel_group = as.character(panel_group), variant_id)]
}

#' Put each gene's own best-variant effect and its decoy best-variant effect
#' on its overview row. Own effects come from the DE high, DE low and
#' positive-control panels; decoy effects from the negative-control panels.
attach_effects <- function(rows, effects) {
  e <- data.table::as.data.table(effects)
  own <- e[panel_group != PANEL_GROUPS[4],
           .(Gene_name, own_slope = slope, own_se = se, own_p = p, own_variant = variant_id)]
  dec <- e[panel_group == PANEL_GROUPS[4],
           .(Gene_name, decoy_slope = slope, decoy_se = se, decoy_p = p)]
  out <- merge(data.table::as.data.table(rows), unique(own, by = "Gene_name"),
               by = "Gene_name", all.x = TRUE, sort = FALSE)
  out <- merge(out, unique(dec, by = "Gene_name"), by = "Gene_name", all.x = TRUE, sort = FALSE)
  out[]
}
