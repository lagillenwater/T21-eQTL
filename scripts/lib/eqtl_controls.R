# eqtl_controls.R
#
# Standalone positive and negative controls for the within-T21 gene-level
# cis-eQTL test, plus the shared runner that the observed, negative and
# positive gene sets all pass through so the three are tested identically.
#
#   gene_tss                 TSS per gene from GTEx allpairs rows
#   select_positive_controls strong GTEx eGenes outside the deviating set
#   assign_decoys            unlinked (distant) variant-set pairing
#   gene_level_tests         per-gene best-variant permutation test + BH
#
# Depends on fit_variants(), perm_min_p() and gene_level_p() from eqtl_fit.R.

#' TSS per gene from allpairs rows: GTEx reports tss_distance = POS - TSS.
#' @param ap data.table with ensembl_stable, POS, tss_distance
#' @return data.table(ensembl_stable, tss); the median over a gene's rows
gene_tss <- function(ap) {
  stopifnot(all(c("ensembl_stable", "POS", "tss_distance") %in% names(ap)))
  ap[, .(tss = stats::median(as.numeric(POS) - as.numeric(tss_distance))),
     by = ensembl_stable]
}

#' Positive-control genes: the eligible candidates with the smallest GTEx
#' nominal p, excluding the ids in `exclude` (the deviating set).
#' @param candidates data.table with ensembl_stable and Gene_name (already
#'   filtered to the expressed, non-repeat, non-deviating lane)
#' @param gtex_min_p data.table(ensembl_stable, gtex_min_p); genes absent
#'   from it have no GTEx coverage and are skipped
#' @return the top-n candidate rows with gtex_min_p, strongest first
select_positive_controls <- function(candidates, gtex_min_p, exclude, n) {
  stopifnot(all(c("ensembl_stable", "Gene_name") %in% names(candidates)),
            all(c("ensembl_stable", "gtex_min_p") %in% names(gtex_min_p)))
  pool <- merge(candidates[!ensembl_stable %in% exclude],
                gtex_min_p[, .(ensembl_stable, gtex_min_p)],
                by = "ensembl_stable")
  data.table::setorder(pool, gtex_min_p)
  utils::head(pool, n)
}

#' Decoy pairing: for each gene, another tested gene whose TSS is at least
#' min_distance away (so its cis variant set cannot be in LD with the gene's
#' own locus), choosing the candidate whose variant count is closest to the
#' gene's own so the decoy test carries the same multiplicity; ties go to the
#' farthest candidate. NA when no gene is far enough.
#' @param tss        named numeric vector of TSS positions (names = gene names)
#' @param n_variants named numeric vector of tested-variant counts, same names
#' @return data.table(Gene_name, decoy_gene, distance)
assign_decoys <- function(tss, n_variants, min_distance) {
  stopifnot(is.numeric(tss), !is.null(names(tss)), length(tss) >= 1,
            is.numeric(n_variants), setequal(names(n_variants), names(tss)))
  nv <- n_variants[names(tss)]
  D  <- abs(outer(tss, tss, "-"))
  eligible <- D >= min_distance
  diag(eligible) <- FALSE
  # key: variant-count gap first, then distance (larger wins) as a fractional
  # tiebreak that can never outweigh a one-variant difference in the gap.
  key <- abs(outer(nv, nv, "-")) - D / (max(D) + 1)
  key[!eligible] <- Inf
  j  <- max.col(-key, ties.method = "first")
  ok <- rowSums(eligible) > 0
  data.table::data.table(
    Gene_name  = names(tss),
    decoy_gene = ifelse(ok, names(tss)[j], NA_character_),
    distance   = ifelse(ok, D[cbind(seq_along(tss), j)], NA_real_))
}

#' Gene-level permutation test for a set of genes.
#'
#' For each gene: fit every variant in variants_of[[gene]] (columns of G_all)
#' against expr_of(gene), take the smallest p, and compare it with the
#' minimum p across the same variants under n_perm permutations of
#' expression (perm_min_p). BH across the set; detected at q < fdr.
#' The per-gene seed is seed_base + position in `genes`, so a set is
#' reproducible and independent of the other sets.
#' @return data.table(Gene_name, n_variants, min_p_obs, best_variant,
#'   p_gene_perm, q_gene_bh, detected)
gene_level_tests <- function(genes, variants_of, G_all, expr_of, n_perm, seed_base, fdr) {
  stopifnot(is.character(genes), is.list(variants_of), is.matrix(G_all),
            is.function(expr_of))
  na_row <- function(g, k) data.table::data.table(
    Gene_name = g, n_variants = k, min_p_obs = NA_real_,
    best_variant = NA_character_, p_gene_perm = NA_real_)
  res <- data.table::rbindlist(lapply(seq_along(genes), function(i) {
    g    <- genes[i]
    cols <- intersect(variants_of[[g]], colnames(G_all))
    if (length(cols) == 0) return(na_row(g, 0L))
    e <- expr_of(g)
    if (all(is.na(e))) return(na_row(g, length(cols)))
    G    <- G_all[, cols, drop = FALSE]
    fits <- fit_variants(G, e)
    # every variant monomorphic/untestable: which.min would return integer(0)
    if (!any(is.finite(fits$p))) return(na_row(g, length(cols)))
    best <- which.min(fits$p)
    data.table::data.table(
      Gene_name = g, n_variants = length(cols),
      min_p_obs = fits$p[best], best_variant = cols[best],
      p_gene_perm = gene_level_p(fits$p[best],
                                 perm_min_p(G, e, n_perm = n_perm,
                                            seed = as.integer(seed_base) + i)))
  }))
  res[, q_gene_bh := stats::p.adjust(p_gene_perm, "BH")]
  res[, detected := !is.na(q_gene_bh) & q_gene_bh < fdr]
  res[]
}
