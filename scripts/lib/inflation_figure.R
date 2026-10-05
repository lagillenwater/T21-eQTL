# inflation_figure.R
#
# The per-run figure of script 17, built from its saved tables, so it can be
# redrawn (Rscript scripts/17_distal_inflation.R --run <name> --figure-only)
# without rerunning the permutations.
#
#   A. QQ of the distal tests on long-arm variants, chr21 assessable genes and
#      the chr22 reference, shaded permutation range, zoomed to expected 0-4
#      with the full range inset (as in the S3 supplement figure).
#   B. Per-gene lambda, chr21 and chr22 columns, deviating genes in their own
#      colour; filled = outside the gene's own permutation 95% range.
#   C. lambda by variant set, adjusted and rank inverse-normal expression.
#   D. Runs with covariates: each covariate against chr21 genotype; filled =
#      outside its own permutation 95% range.

INFLATION_COL_CHR <- c(chr21 = "#B2182B", chr22 = "#2166AC")
INFLATION_COL_DEV <- "#1B9E77"   # deviating genes; distinct from the run colours of S3

#' Per-gene lambda, one column per chromosome, deviating genes in their own
#' colour; filled = outside the gene's own permutation 95% range (each gene
#' has its own range, so no shared band is drawn). Shared by script 17 and S3.
#' @param gene_out script 17's distal_inflation_genes.csv
#' @param col_labels column label per chromosome, c(chr21 =, chr22 =)
#' @param label_outside also label chr21 genes outside their own range
per_gene_lambda_panel <- function(gene_out, col_labels, title, subtitle, seed = 1717L, label_outside = TRUE) {
  g <- data.table::copy(data.table::as.data.table(gene_out))
  g[, col_f := factor(col_labels[chromosome], levels = col_labels)]
  g[, deviating := gene_set == "Deviating"]
  # Label the deviating genes, the five highest lambdas overall and, if
  # asked, the chr21 genes outside their own range.
  g[, label := ifelse(deviating | (label_outside & above_range & chromosome == "chr21") |
                        rank(-lambda, ties.method = "first") <= 5, Gene_name, "")]
  data.table::setorder(g, deviating)                 # deviating genes drawn on top
  jit <- position_jitter(width = 0.3, height = 0, seed = seed)
  ggplot(g, aes(col_f, lambda)) +
    geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_point(aes(colour = deviating, shape = above_range, size = deviating), position = jit, alpha = 0.9,
               stroke = 0.6) +
    ggrepel::geom_text_repel(aes(label = label), position = jit, size = 2.3, colour = "grey15",
                             max.overlaps = 40, segment.size = 0.2, min.segment.length = 0) +
    scale_colour_manual(values = c(`TRUE` = INFLATION_COL_DEV, `FALSE` = "grey55"),
                        labels = c(`TRUE` = "Deviating gene", `FALSE` = "Other gene"), name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = "Outside its own permutation range",
                                  `FALSE` = "Inside its own permutation range"), name = NULL) +
    scale_size_manual(values = c(`TRUE` = 2.3, `FALSE` = 1.2), guide = "none") +
    labs(x = NULL, y = "Per-gene lambda", title = title, subtitle = subtitle) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.direction = "vertical")
}

#' @param summary_tbl,qq_tbl,gene_out,cov_out the tables script 17 writes
#'   (cov_out NULL for a run without covariates)
#' @param control_set the chr22 set name
#' @return list(fig, height)
distal_inflation_figure <- function(summary_tbl, qq_tbl, gene_out, cov_out, control_set,
                                    min_dist, seed = 1717L, zoom_x = 4) {
  # Panel A matches the S3 supplement figure: long-arm variants (the short
  # arm holds GRCh38's false duplications), shaded permutation range, zoomed
  # to expected 0 to zoom_x where any inflation shows, full range inset.
  key <- c("gene_set", "variant_set")
  qq_tbl <- merge(data.table::as.data.table(qq_tbl)[variant_set == "q-arm only"],
                  summary_tbl[variant_set == "q-arm only",
                              .(gene_set, variant_set, n_genes, lambda, perm_lo_l = perm_lo, perm_hi_l = perm_hi)],
                  by = key)
  qq_tbl[, curve := sprintf("%s (%d genes): lambda %.3f (permutation 95%% %.3f-%.3f)",
                            gene_set, n_genes, lambda, perm_lo_l, perm_hi_l)]
  curves <- unique(qq_tbl[, .(curve, gene_set)])
  col_curve <- setNames(INFLATION_COL_CHR[ifelse(curves$gene_set == control_set, "chr22", "chr21")], curves$curve)
  x_max <- max(qq_tbl$expected) * 1.02
  y_max <- max(qq_tbl$perm_hi, qq_tbl$observed) * 1.02
  zoom_y <- max(qq_tbl[expected <= zoom_x, c(observed, perm_hi)]) * 1.02
  qq_layers <- function(d) list(
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4),
    geom_ribbon(data = d, aes(ymin = perm_lo, ymax = perm_hi, fill = curve), colour = NA, alpha = 0.18),
    geom_line(data = d, linewidth = 0.6),
    scale_colour_manual(values = col_curve, name = NULL),
    scale_fill_manual(values = col_curve, guide = "none"))
  zoomed <- qq_tbl[expected <= zoom_x * 1.05]
  main <- ggplot(zoomed, aes(expected, observed, colour = curve)) +
    qq_layers(zoomed) +
    coord_cartesian(xlim = c(0, zoom_x), ylim = c(0, zoom_y), expand = FALSE) +
    labs(x = expression(Expected ~ -log[10] ~ italic(p)), y = expression(Observed ~ -log[10] ~ italic(p)),
         title = "A  Distal genotype-expression tests, T21 subjects",
         subtitle = sprintf("Common chr21 long-arm variants >= %.0f Mb from each TSS; shaded: permutation 95%% range",
                            min_dist / 1e6)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.minor = element_blank(),
          aspect.ratio = 1)
  inset <- ggplot(qq_tbl, aes(expected, observed, colour = curve)) +
    qq_layers(qq_tbl) +
    annotate("rect", xmin = 0, xmax = zoom_x, ymin = 0, ymax = zoom_y, fill = NA, colour = "grey30",
             linewidth = 0.3) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE) +
    guides(colour = "none") +
    labs(x = NULL, y = NULL, title = "full range") +
    theme_bw(base_size = 6) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 6),
          plot.background = element_rect(colour = "grey70"))
  pA <- main + inset_element(inset, left = 0.55, bottom = 0.03, right = 0.98, top = 0.45)

  pB <- per_gene_lambda_panel(gene_out, c(chr21 = "chr21 assessable", chr22 = control_set),
                              "B  Per-gene lambda",
                              "Long-arm variants; filled: outside the gene's own permutation 95% range", seed)

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
