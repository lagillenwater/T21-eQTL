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
