# inflation_figure.R
#
# The per-run figure of script 17, built from its saved tables, so it can be
# redrawn (Rscript scripts/17_distal_inflation.R --run <name> --figure-only)
# without rerunning the permutations.
#
#   A. QQ of the distal tests, all variants (solid, shaded permutation range)
#      and long-arm variants (dashed), chr21 assessable genes and the chr22
#      reference.
#   B. Per-gene lambda by gene set; filled = outside the gene's own
#      permutation 95% range (each gene has its own range, so no shared band
#      is drawn).
#   C. lambda by variant set, adjusted and rank inverse-normal expression.
#   D. Runs with covariates: each covariate against chr21 genotype; filled =
#      outside its own permutation 95% range.

INFLATION_COL_CHR <- c(chr21 = "#B2182B", chr22 = "#2166AC")

#' @param summary_tbl,qq_tbl,gene_out,cov_out the tables script 17 writes
#'   (cov_out NULL for a run without covariates)
#' @param groups gene-set order for panel B; control_set the chr22 set name
#' @return list(fig, height)
distal_inflation_figure <- function(summary_tbl, qq_tbl, gene_out, cov_out, groups, control_set,
                                    min_dist, seed = 1717L) {
  key <- c("gene_set", "variant_set")
  qq_tbl <- merge(data.table::as.data.table(qq_tbl),
                  summary_tbl[, .(gene_set, variant_set, n_genes, lambda, perm_lo_l = perm_lo,
                                  perm_hi_l = perm_hi)], by = key)
  qq_tbl[, curve := sprintf("%s (%d genes), %s variants: lambda %.3f (permutation 95%% %.3f-%.3f)",
                            gene_set, n_genes, sub(" only", "", variant_set), lambda, perm_lo_l, perm_hi_l)]
  curves <- unique(qq_tbl[, .(curve, gene_set, variant_set)])
  col_curve <- setNames(INFLATION_COL_CHR[ifelse(curves$gene_set == control_set, "chr22", "chr21")], curves$curve)
  lty_curve <- setNames(ifelse(curves$variant_set == "all", "solid", "22"), curves$curve)
  # The axis is set by the long-arm curves and the envelopes; the short-arm
  # false-duplication tail of the all-variant chr21 curve runs off the top.
  y_cap <- ceiling(max(qq_tbl[variant_set == "q-arm only" | gene_set == control_set, observed],
                       qq_tbl$perm_hi, qq_tbl$expected)) + 1
  off_top <- qq_tbl[observed > y_cap, max(observed)]
  pA <- ggplot(qq_tbl, aes(expected, observed, colour = curve, linetype = curve)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_ribbon(data = qq_tbl[variant_set == "all"], aes(ymin = perm_lo, ymax = perm_hi, fill = curve),
                colour = NA, alpha = 0.18, show.legend = FALSE) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = col_curve, name = NULL) +
    scale_fill_manual(values = col_curve, guide = "none") +
    scale_linetype_manual(values = lty_curve, name = NULL) +
    coord_cartesian(xlim = c(0, y_cap), ylim = c(0, y_cap), expand = FALSE) +
    labs(x = expression(Expected ~ -log[10] ~ italic(p)), y = expression(Observed ~ -log[10] ~ italic(p)),
         title = "A  Distal genotype-expression tests, T21 subjects",
         subtitle = sprintf("Common chr21 variants >= %.0f Mb from each TSS; shaded: permutation 95%% range",
                            min_dist / 1e6)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.key.width = unit(1.4, "lines"),
          panel.grid.minor = element_blank(), aspect.ratio = 1)
  if (length(off_top) && is.finite(off_top)) {
    pA <- pA + annotate("text", x = 0.3, y = y_cap - 0.4, hjust = 0, vjust = 1, size = 2.5, colour = "grey25",
                        label = "chr21, all variants: leaves the panel here\n(short-arm false-duplication variants)")
  }

  g <- data.table::copy(data.table::as.data.table(gene_out))
  g[, set_f := factor(gene_set, levels = groups[groups %in% gene_set])]
  # Label the deviating genes, the chr21 genes outside their own range, and
  # the five highest lambdas overall.
  g[, label := ifelse(gene_set == "Deviating" | (above_range & chromosome == "chr21") |
                        rank(-lambda, ties.method = "first") <= 5, Gene_name, "")]
  pB <- ggplot(g, aes(set_f, lambda, colour = chromosome)) +
    geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_point(aes(shape = above_range), position = position_jitter(width = 0.25, height = 0, seed = seed),
               size = 1.3, alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.3, colour = "grey15", max.overlaps = 40,
                             position = position_jitter(width = 0.25, height = 0, seed = seed),
                             segment.size = 0.2, min.segment.length = 0) +
    scale_colour_manual(values = INFLATION_COL_CHR, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    scale_x_discrete(labels = function(x) sub(", ", ",\n", x)) +
    labs(x = NULL, y = "Per-gene lambda (long-arm variants)",
         title = "B  Per-gene lambda",
         subtitle = "Filled: outside the gene's own permutation 95% range") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())

  var_tbl <- summary_tbl[gene_set %in% c("chr21 assessable", control_set)]
  var_tbl[, var_f := factor(variant_set, levels = unique(summary_tbl$variant_set))]
  var_tbl[, chromosome := ifelse(gene_set == control_set, "chr22", "chr21")]
  dodge <- position_dodge(width = 0.5)
  pC <- ggplot(var_tbl, aes(var_f, lambda, colour = chromosome, group = chromosome)) +
    geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(ymin = perm_lo, ymax = perm_hi), position = dodge, linewidth = 2.2, colour = "grey80") +
    geom_point(position = dodge, size = 2) +
    geom_point(aes(y = lambda_inverse_normal), position = dodge, size = 2, shape = 1) +
    scale_colour_manual(values = INFLATION_COL_CHR, labels = c(chr21 = "chr21 assessable", chr22 = control_set),
                        name = NULL) +
    labs(x = "Variants (T21 minor allele frequency bins use both arms)", y = "lambda",
         title = "C  lambda by variant set",
         subtitle = "Filled: adjusted expression; open: rank inverse-normal; grey: permutation 95% range") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())

  if (is.null(cov_out)) return(list(fig = pA + (pB / pC) + plot_layout(widths = c(1.1, 1)), height = 7.5))
  cv <- data.table::copy(data.table::as.data.table(cov_out))
  cv[, covariate_f := factor(display_covariate(covariate), levels = rev(display_covariate(covariate)))]
  pD <- ggplot(cv, aes(lambda, covariate_f)) +
    geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(xmin = perm_lo, xmax = perm_hi), linewidth = 2.2, colour = "grey80") +
    geom_point(aes(shape = above_range), size = 1.8, colour = "grey15") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(x = "lambda (covariate ~ each common chr21 long-arm variant)", y = NULL,
         title = "D  Run covariates against chr21 genotype",
         subtitle = "Filled: outside its own permutation 95% range (grey)") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
  list(fig = (pA | pB) / (pC | pD), height = 12)
}
