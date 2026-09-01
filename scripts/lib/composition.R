# composition.R
#
# Composition control for chr21 deviating genes.
#
# A chr21 gene can appear to deviate from dosage expectation because the
# blood cell type that expresses it changed in abundance in T21, not because
# of any regulatory effect on the gene itself. Composition shifts move whole
# co-expression programs, not single genes, so the diagnostic is: does this
# gene's co-expression neighborhood (its 20 most correlated non-chr21
# partners) shift with it in T21 vs Control?
#
# Partners are found in CONTROLS ONLY so the karyotype effect cannot leak
# into the correlation used to define the neighborhood - correlating across
# both karyotypes would let a shared T21-vs-Control shift inflate the
# correlation of any two genes that both respond to trisomy, which is
# exactly the confound this control is meant to catch.
#
# The null is the median log2FC of 2000 random 20-gene sets drawn from the
# same background (expressed, non-chr21 genes) - "how much would a random
# gene set of this size appear to shift, just from background noise/other
# programs".

#' Median log2FC of a gene's top-n co-expressed non-chr21 partners.
#'
#' @param gene        gene name, must be a row name of L_ctrl and present in lfc
#' @param L_ctrl      genes x controls log2-CPM matrix of expressed, non-chr21
#'                    genes (the partner pool)
#' @param gene_ctrl   the target gene's log2-CPM across the same controls
#'                    (same column order as L_ctrl)
#' @param lfc         named vector of genome-wide T21-vs-Control log2FC
#' @param n_partners  number of top-correlated partners to use
#' @return scalar median log2FC of the top n_partners genes
partner_shift <- function(gene, L_ctrl, gene_ctrl, lfc, n_partners = 20) {
  stopifnot(is.matrix(L_ctrl), is.numeric(gene_ctrl), is.numeric(lfc))
  r <- cor(gene_ctrl, t(L_ctrl))[1, ]
  top <- names(sort(r, decreasing = TRUE))[seq_len(n_partners)]
  median(lfc[top], na.rm = TRUE)
}

#' Null distribution for partner_shift: median log2FC of random gene sets.
#'
#' Vectorised - draws all n_draw * n_partners indices at once and takes
#' column medians, rather than looping over draws.
#'
#' @param lfc_bg    background vector of log2FC values to sample from
#' @param n_partners  genes per random set
#' @param n_draw    number of random sets to draw
#' @param seed      RNG seed, for reproducibility
#' @return numeric vector of length n_draw, the median log2FC of each set
partner_null <- function(lfc_bg, n_partners = 20, n_draw = 2000, seed = 1) {
  stopifnot(is.numeric(lfc_bg))
  set.seed(seed)
  idx <- matrix(sample.int(length(lfc_bg), n_partners * n_draw, replace = TRUE),
                nrow = n_partners, ncol = n_draw)
  draws <- matrix(lfc_bg[idx], nrow = n_partners, ncol = n_draw)
  apply(draws, 2, median, na.rm = TRUE)
}

#' Classify a gene's deviation as program-driven, mixed, or gene-specific.
#'
#' The empirical p is one-sided in the gene's own deviation direction: how
#' often does a random gene set shift at least as far as the observed
#' partner shift, in the direction the gene itself moved.
#'
#' @param gene_lfc     the gene's own T21-vs-Control log2FC
#' @param partner_lfc  median log2FC of the gene's co-expressed partners
#' @param p            empirical p-value (see partner_null / the p computed
#'                      by the caller as mean(sign(gene_lfc) * null >=
#'                      sign(gene_lfc) * partner_lfc))
#' @return one of "PROGRAM", "MIXED", "GENE-SPECIFIC"
composition_verdict <- function(gene_lfc, partner_lfc, p) {
  stopifnot(is.numeric(gene_lfc), is.numeric(partner_lfc), is.numeric(p))
  share <- partner_lfc / gene_lfc
  if (p < 0.05 && share >= 0.5) {
    "PROGRAM"
  } else if (p < 0.05) {
    "MIXED"
  } else {
    "GENE-SPECIFIC"
  }
}
