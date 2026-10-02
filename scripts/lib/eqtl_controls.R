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

# =============================================================================
# Controls v2: matched, one-per-locus positive controls; several decoy sets
# per gene with a MAF floor (docs/decisions.md, controls v2).
# =============================================================================

#' Positive-control genes matched to the deviating genes.
#'
#' One control per deviating gene, chosen from GTEx whole-blood eGenes among
#' the expressed, non-repeat, non-deviating chr21 genes: the candidate nearest
#' the deviating gene in standardised (log2 |aFC|, log10 baseMean) space, so
#' the controls carry effect sizes and expression levels like the genes they
#' stand in for rather than the strongest eQTLs on the chromosome. Greedy in
#' order of decreasing target |aFC|; a candidate is used once, and its TSS
#' must be at least `min_separation` from every control already chosen, so
#' two controls never share a locus.
#' @param targets    data.table(Gene_name, abs_afc, baseMean) of the deviating
#'   genes that have a GTEx aFC; rows with NA are skipped
#' @param candidates data.table(ensembl_stable, Gene_name, abs_afc, baseMean,
#'   tss), already filtered to eligible eGenes with exclusions removed
#' @return data.table(target_gene, target_abs_afc, target_baseMean,
#'   Gene_name, ensembl_stable, abs_afc, baseMean, tss, match_distance)
match_positive_controls <- function(targets, candidates, min_separation) {
  stopifnot(all(c("Gene_name", "abs_afc", "baseMean") %in% names(targets)),
            all(c("ensembl_stable", "Gene_name", "abs_afc", "baseMean", "tss") %in% names(candidates)))
  tg <- data.table::as.data.table(targets)[!is.na(abs_afc) & !is.na(baseMean) & abs_afc > 0 & baseMean > 0]
  cd <- data.table::as.data.table(candidates)[!is.na(abs_afc) & !is.na(baseMean) & !is.na(tss) &
                                                abs_afc > 0 & baseMean > 0]
  data.table::setorder(tg, -abs_afc)
  z <- function(x, ref) { s <- stats::sd(ref); (x - mean(ref)) / (if (is.na(s) || s == 0) 1 else s) }
  ca <- log2(cd$abs_afc); cb <- log10(cd$baseMean)
  cz <- cbind(z(ca, ca), z(cb, cb))
  tz <- cbind(z(log2(tg$abs_afc), ca), z(log10(tg$baseMean), cb))
  # distance matrix targets x candidates, computed once
  D <- sqrt(outer(tz[, 1], cz[, 1], "-")^2 + outer(tz[, 2], cz[, 2], "-")^2)
  avail <- rep(TRUE, nrow(cd)); chosen <- integer(0)
  # greedy over targets: each pick removes one candidate and its locus
  for (i in seq_len(nrow(tg))) {
    d <- D[i, ]; d[!avail] <- Inf
    if (all(is.infinite(d))) { chosen <- c(chosen, NA_integer_); next }
    j <- which.min(d)
    chosen <- c(chosen, j)
    avail[abs(cd$tss - cd$tss[j]) < min_separation] <- FALSE
  }
  out <- data.table::data.table(target_gene = tg$Gene_name, target_abs_afc = tg$abs_afc,
                                target_baseMean = tg$baseMean)
  pick <- cd[ifelse(is.na(chosen), NA_integer_, chosen),
             .(Gene_name, ensembl_stable, abs_afc, baseMean, tss)]
  out <- cbind(out, pick)
  out[, match_distance := D[cbind(seq_len(nrow(out)), chosen)]]
  out[!is.na(Gene_name)]
}

#' Several decoy variant sets per gene.
#'
#' Like assign_decoys, but returns the `k` best donors per gene (closest
#' variant count first, farther wins ties) among those at least
#' `min_distance` away, ranked 1..k. Donors may include genes that are not
#' themselves being tested against decoys (e.g. positive controls), which is
#' what `donor_tss` / `donor_n_variants` allow.
#' @param tss,n_variants     named vectors for the genes to test
#' @param donor_tss,donor_n_variants named vectors for the donor pool
#' @return data.table(Gene_name, decoy_gene, decoy_rank, distance, n_variants)
assign_decoy_sets <- function(tss, n_variants, donor_tss, donor_n_variants, min_distance, k) {
  stopifnot(is.numeric(tss), !is.null(names(tss)), setequal(names(n_variants), names(tss)),
            is.numeric(donor_tss), !is.null(names(donor_tss)),
            setequal(names(donor_n_variants), names(donor_tss)), k >= 1)
  nv  <- n_variants[names(tss)]; dnv <- donor_n_variants[names(donor_tss)]
  D   <- abs(outer(tss, donor_tss, "-"))
  rank_key <- abs(outer(nv, dnv, "-")) - D / (max(D) + 1)
  long <- data.table::data.table(Gene_name = rep(names(tss), times = length(donor_tss)),
                                 decoy_gene = rep(names(donor_tss), each = length(tss)),
                                 distance = as.vector(D), rank_key = as.vector(rank_key),
                                 n_variants = rep(unname(dnv), each = length(tss)))
  long <- long[Gene_name != decoy_gene & distance >= min_distance]
  data.table::setorder(long, Gene_name, rank_key)
  long[, decoy_rank := seq_len(.N), by = Gene_name]
  long[decoy_rank <= k, .(Gene_name, decoy_gene, decoy_rank, distance, n_variants)]
}
