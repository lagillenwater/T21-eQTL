# inflation.R
#
# Genomic inflation of distal genotype-expression tests: does the T21
# genotype-expression data carry structure (ancestry, relatedness, batch)
# that would inflate the within-T21 eQTL tests?
#
# Every expression gene is regressed on every common chr21 variant far enough
# from its TSS that the variant cannot tag the gene's own cis locus: one
# ordinary least-squares fit per (gene, variant) pair, the same fit as
# fit_variants() in eqtl_fit.R. Without structure those p-values are uniform
# and the genomic inflation factor lambda (median chi-square over its null
# median) is 1.
#
# A full scan is tens of millions of pairs, so nothing per pair is kept. Each
# |r| is binned on an exact -log10 p grid (bin edges precomputed in |r|, so no
# pt() call per pair) and only the bin counts are accumulated, per gene group
# and MAF bin; lambda and the QQ curve are read off the counts. Per-gene
# lambda comes from the per-gene median |r|, which maps to the median p
# exactly because p is monotone in |r|.
#
# The reference for lambda is a permutation of subjects applied jointly to
# every gene: it keeps co-expression between genes and LD between variants,
# so its spread is the sampling spread of lambda under these correlated
# tests, which no closed form gives.

#' Genomic inflation factor from p-values (1 df chi-square scale).
lambda_from_p <- function(p) {
  qchisq(median(p, na.rm = TRUE), 1, lower.tail = FALSE) / qchisq(0.5, 1)
}

#' p-value of the OLS slope test from a correlation, with df residual degrees.
p_from_abs_r <- function(r, df) {
  t <- r * sqrt(df / pmax(1 - r^2, .Machine$double.eps))
  2 * pt(-t, df)
}

#' Lambda from a median |r| (the median p is the p of the median |r|).
lambda_from_abs_r <- function(r_med, df) {
  qchisq(p_from_abs_r(r_med, df), 1, lower.tail = FALSE) / qchisq(0.5, 1)
}

# project_out() and covariate_rank() live in eqtl_fit.R, shared with the
# eQTL fits; source it first.

#' Centre each column, set missing entries to the column mean (0 after
#' centring) and scale to unit length, so crossprod() of two such matrices is
#' the matrix of Pearson correlations. Constant columns come back all zero.
standardise_cols <- function(M, covariates = NULL) {
  mu <- colMeans(M, na.rm = TRUE)
  M  <- sweep(M, 2, mu, "-")
  M[is.na(M)] <- 0
  if (!is.null(covariates)) M <- project_out(M, covariates)
  len <- sqrt(colSums(M^2))
  len[len == 0] <- Inf
  M / rep(len, each = nrow(M))
}

#' Rank-based inverse normal transform of each column (ties averaged).
inverse_normal_cols <- function(M) {
  R <- apply(M, 2, rank, ties.method = "average", na.last = "keep")
  n <- colSums(!is.na(M))
  qnorm((R - 0.5) / rep(n, each = nrow(M)))
}

#' Bin grid on -log10 p and the matching |r| edges.
#'
#' Bin k holds -log10 p in [x[k], x[k + 1]); the last bin is open above.
#' @param df residual degrees of freedom of each test (n - 2 - n_covariates)
#' @return list(x = lower -log10 p edge per bin, r = matching |r| edge, df)
log10p_bins <- function(df, width = 1e-4, max_x = 50) {
  x <- seq(0, max_x, by = width)
  t <- qt(10^(-x) / 2, df, lower.tail = FALSE)
  r <- t / sqrt(df + t^2)
  r[1] <- 0
  list(x = x, r = r, df = df)
}

#' Per-bin counts: an array from distal_scan() is summed over its other
#' dimensions (MAF bins, gene groups) first.
bin_counts <- function(counts) {
  if (is.null(dim(counts))) as.numeric(counts) else rowSums(counts)
}

#' Median -log10 p and lambda from bin counts.
hist_lambda <- function(counts, bins) {
  counts <- bin_counts(counts)
  n <- sum(counts)
  if (n == 0) return(NA_real_)
  k <- which(cumsum(counts) >= n / 2)[1]
  # Bin midpoint; the half-width (5e-5 on -log10 p) moves lambda by < 0.1%.
  x <- bins$x[k] + diff(bins$x[1:2]) / 2
  qchisq(10^(-x), 1, lower.tail = FALSE) / qchisq(0.5, 1)
}

#' QQ curve from bin counts at log-spaced ranks (smallest p first).
#' @return data.table(rank, expected, observed) on the -log10 p scale
hist_qq <- function(counts, bins, n_points = 400) {
  counts <- bin_counts(counts)
  n <- sum(counts)
  if (n == 0) return(data.table::data.table(rank = numeric(0), expected = numeric(0),
                                            observed = numeric(0)))
  ranks <- unique(round(10^seq(0, log10(n), length.out = n_points)))
  from_top <- rev(cumsum(rev(counts)))   # tests at or above bin k
  # rank-th smallest p lies in the highest bin whose count-from-top reaches rank
  k <- findInterval(-ranks, -from_top)
  data.table::data.table(rank = ranks, expected = -log10((ranks - 0.5) / n),
                         observed = bins$x[pmax(k, 1)])
}

#' One pass of the distal scan.
#'
#' @param Ys n x g standardised expression (standardise_cols)
#' @param Gs n x m standardised genotypes, same subjects and covariates
#' @param gene_tss TSS per gene (NA = different chromosome: every variant is distal)
#' @param pos position per variant
#' @param gene_group integer group per gene (1..n_groups)
#' @param variant_bin integer MAF bin per variant (1..n_bins_maf)
#' @param min_dist minimum |pos - TSS| for a pair to count as distal
#' @param bins log10p_bins() grid
#' @param hit_x record pairs at -log10 p >= hit_x (NULL: record none)
#' @param stat_bins variant categories the per-gene median, maximum and test
#'   count use (NULL: all); the bin counts and hits always cover every category
#' @return list(counts = array [p bin, maf bin, group], gene_median_r,
#'   gene_n, gene_max_r, hits = data.table(gene, variant, r) or NULL)
distal_scan <- function(Ys, Gs, gene_tss, pos, gene_group, variant_bin, min_dist,
                        bins, n_groups = max(gene_group), n_maf = max(variant_bin),
                        hit_x = NULL, block = 50L, stat_bins = NULL) {
  stopifnot(nrow(Ys) == nrow(Gs), length(gene_tss) == ncol(Ys), length(pos) == ncol(Gs))
  nb  <- length(bins$x)
  g   <- ncol(Ys)
  counts <- numeric(nb * n_maf * n_groups)
  med <- n_test <- max_r <- rep(NA_real_, g)
  hits <- list()
  hit_r <- if (is.null(hit_x)) Inf else bins$r[findInterval(hit_x, bins$x)]
  m <- ncol(Gs)
  in_stat <- if (is.null(stat_bins)) rep(TRUE, m) else variant_bin %in% stat_bins
  for (idx in split(seq_len(g), ceiling(seq_len(g) / block))) {
    A <- abs(crossprod(Ys[, idx, drop = FALSE], Gs))       # |r|, genes x variants
    tss <- gene_tss[idx]
    keep <- if (all(is.na(tss))) {
      matrix(TRUE, length(idx), m)
    } else {
      dist <- abs(outer(ifelse(is.na(tss), -Inf, tss), pos, "-"))
      dist >= min_dist
    }
    a  <- A[keep]
    vb <- variant_bin[col(A)[keep]]
    gb <- gene_group[idx][row(A)[keep]]
    cell <- findInterval(a, bins$r) + nb * (vb - 1L) + nb * n_maf * (gb - 1L)
    counts <- counts + tabulate(cell, nbins = length(counts))
    A[!keep] <- NA
    if (is.finite(hit_r)) {
      w <- which(A >= hit_r, arr.ind = TRUE)
      if (nrow(w)) hits[[length(hits) + 1L]] <-
        data.table::data.table(gene = idx[w[, 1]], variant = w[, 2], r = A[w])
    }
    A[, !in_stat] <- NA
    med[idx]    <- matrixStats::rowMedians(A, na.rm = TRUE)
    max_r[idx]  <- matrixStats::rowMaxs(A, na.rm = TRUE)
    n_test[idx] <- rowSums(keep[, in_stat, drop = FALSE])
  }
  list(counts = array(counts, c(nb, n_maf, n_groups)),
       gene_median_r = med, gene_n = n_test, gene_max_r = max_r,
       hits = if (length(hits)) data.table::rbindlist(hits) else NULL)
}

#' Variant QC for a T21 dosage matrix (subjects x variants).
#' @return data.table(call_rate, alt_af, maf) per variant
dosage_stats <- function(G, ploidy = 3) {
  called <- colSums(!is.na(G))
  af <- colSums(G, na.rm = TRUE) / (ploidy * pmax(called, 1))
  data.table::data.table(call_rate = called / nrow(G), alt_af = af, maf = pmin(af, 1 - af))
}

#' Collapse a variants x samples dosage table to a subjects x variants matrix.
#'
#' Sample columns are mapped to subjects with subject_id_from_labid()
#' (cohort.R); a subject with several sample columns gets the mean of its
#' called dosages, as in scripts 03 and 16.
#' @param dt data.table with POS, REF, ALT and one column per sample
#' @param subjects subject ids wanted, in row order of the result
dosage_matrix <- function(dt, subjects) {
  sample_cols <- setdiff(names(dt), c("POS", "REF", "ALT"))
  subj <- subject_id_from_labid(sample_cols)
  use  <- subj %in% subjects
  M  <- t(as.matrix(dt[, sample_cols[use], with = FALSE]))   # samples x variants
  grp <- factor(subj[use], levels = subjects)
  num <- rowsum(ifelse(is.na(M), 0, M), grp, reorder = TRUE)
  den <- rowsum((!is.na(M)) * 1, grp, reorder = TRUE)
  G <- num / den
  G[den == 0] <- NA
  G <- G[match(subjects, rownames(G)), , drop = FALSE]
  rownames(G) <- subjects
  colnames(G) <- paste("chr21", dt$POS, dt$REF, dt$ALT, sep = "_")
  G
}

#' Gene TSS (strand-aware gene start) from a GENCODE GTF, for given chromosomes.
gene_tss_gtf <- function(gtf_path, chromosomes = c("chr21", "chr22")) {
  chr_test <- paste(sprintf('$1=="%s"', chromosomes), collapse = " || ")
  cmd <- sprintf("gzip -dc %s | awk -F'\\t' '$3==\"gene\" && (%s)'", shQuote(gtf_path), chr_test)
  gtf <- data.table::fread(cmd = cmd, header = FALSE, sep = "\t", quote = "")
  data.table::data.table(
    ensembl_stable = sub("\\..*$", "", sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)),
    chromosome = gtf$V1,
    tss = ifelse(gtf$V7 == "+", gtf$V4, gtf$V5))
}
