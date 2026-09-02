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
# The null must be CORRELATION-MATCHED. An earlier version drew independent
# random 20-gene sets, which is anti-conservative: the observed partners are a
# co-expression module and move together, so their median log2FC is far more
# variable than an independent set's (measured on this cohort: SD 0.259 for
# correlation-matched draws vs 0.058 for independent ones). The null therefore
# draws random SEED genes from the same background pool and, for each seed,
# takes ITS OWN top-n_partners correlated genes (controls-only correlation,
# exactly as the real computation does) - so a null draw is a co-expression
# module of the same construction as the observed one, and the comparison is
# "does this gene's module shift more than a typical module", not "more than a
# typical random gene list".

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

#' Correlation-matched null distribution for partner_shift.
#'
#' Draws n_draw random SEED genes from the partner pool and, for each seed,
#' returns the median log2FC of that seed's own top-n_partners co-expressed
#' genes (correlation computed in controls only, as in partner_shift). A null
#' draw is therefore a co-expression module built the same way as the observed
#' one, not an independent random gene list.
#'
#' Vectorised: the full seed-by-pool correlation matrix is computed in a single
#' cor() call, then each row's top partners are taken by row-wise ordering.
#' There is no loop over draws re-calling cor().
#'
#' The seed gene is excluded from its own partner set (its self-correlation is
#' 1), matching the real computation where the chr21 target gene is not a
#' member of the non-chr21 partner pool.
#'
#' @param L_ctrl      genes x controls log2-CPM matrix of the partner pool
#'                    (expressed, non-chr21 genes); rownames are gene names
#' @param lfc         named vector of genome-wide T21-vs-Control log2FC
#' @param n_partners  genes per module
#' @param n_draw      number of seed genes to draw
#' @param seed        RNG seed, for reproducibility
#' @return numeric vector of length n_draw, the median partner log2FC per draw
partner_null <- function(L_ctrl, lfc, n_partners = 20, n_draw = 300, seed = 1) {
  stopifnot(is.matrix(L_ctrl), is.numeric(lfc), !is.null(rownames(L_ctrl)))
  n_pool <- nrow(L_ctrl)
  stopifnot(n_pool > n_partners + 1L)

  set.seed(seed)
  seed_idx <- sample.int(n_pool, n_draw, replace = n_draw > n_pool)

  # n_draw x n_pool correlation matrix, one cor() call.
  R <- suppressWarnings(cor(t(L_ctrl[seed_idx, , drop = FALSE]), t(L_ctrl)))

  # Drop each seed's self-correlation so it cannot be its own top partner.
  R[cbind(seq_len(n_draw), seed_idx)] <- -Inf

  lfc_pool <- lfc[rownames(L_ctrl)]
  apply(R, 1, function(r) {
    top <- order(r, decreasing = TRUE, na.last = TRUE)[seq_len(n_partners)]
    median(lfc_pool[top], na.rm = TRUE)
  })
}

#' One-sided empirical p for an observed partner shift against the null.
#'
#' "How often does a null module shift at least as far as the observed partner
#' set, in the direction the gene itself moved." The (1 + k) / (n_draw + 1)
#' form is the standard permutation p-value: it can never be exactly 0, which
#' would otherwise claim more resolution than n_draw draws can supply.
#'
#' @param gene_lfc     the gene's own T21-vs-Control log2FC (supplies direction)
#' @param partner_lfc  observed median log2FC of the gene's partners
#' @param null         numeric vector from partner_null()
#' @return scalar p in (0, 1]
partner_p <- function(gene_lfc, partner_lfc, null) {
  stopifnot(is.numeric(gene_lfc), is.numeric(partner_lfc), is.numeric(null),
            length(null) > 0)
  s <- sign(gene_lfc)
  (1 + sum(s * null >= s * partner_lfc)) / (length(null) + 1)
}

#' Classify a gene's deviation as program-driven, mixed, or gene-specific.
#'
#' The empirical p is one-sided in the gene's own deviation direction: how
#' often does a random gene set shift at least as far as the observed
#' partner shift, in the direction the gene itself moved.
#'
#' @param gene_lfc     the gene's own T21-vs-Control log2FC
#' @param partner_lfc  median log2FC of the gene's co-expressed partners
#' @param p            empirical p-value from partner_p()
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
