# covariates.R
#
# Covariate handling for the run mechanism: the HTP CyTOF cluster-fraction
# reader, composition fractions, the DESeq2 design formula, the
# expression-artifact adjustment (script 01), vectorized least-squares
# coefficients, and the per-covariate attribution of a karyotype-coefficient
# change (S1). All outputs are descriptive: association, not cause.

#' Read the Synapse CyTOF cluster-percentage table into a samples x clusters
#' matrix.
#'
#' @param path long-format file with columns LabID, Cell_cluster_name, Value
#' @param na_as_zero a missing (LabID, cluster) pair is a population that was
#'   not observed in that sample; TRUE records it as 0 percent
#' @return numeric matrix, rownames = LabID, colnames = cluster names
read_cytof_wide <- function(path, na_as_zero = TRUE) {
  long <- data.table::fread(path)
  stopifnot(all(c("LabID", "Cell_cluster_name", "Value") %in% names(long)))
  if (anyDuplicated(long[, .(LabID, Cell_cluster_name)]))
    stop("CyTOF table has more than one row for a (LabID, cluster) pair", call. = FALSE)
  wide <- data.table::dcast(long, LabID ~ Cell_cluster_name, value.var = "Value")
  X <- as.matrix(wide[, -1])
  rownames(X) <- wide$LabID
  if (na_as_zero) X[is.na(X)] <- 0
  X
}

#' Least-squares coefficients of every row of Y on the design X.
#'
#' @param Y genes x samples
#' @param X samples x p design matrix (include the intercept column yourself)
#' @return genes x p coefficient matrix
coef_rows <- function(Y, X) {
  X <- as.matrix(X)
  stopifnot(is.matrix(Y), nrow(X) == ncol(Y))
  B <- Y %*% t(solve(crossprod(X), t(X)))
  colnames(B) <- colnames(X)
  B
}

#' Composition fractions in percent, reference cluster dropped, rows in the
#' order of `lab_ids`.
composition_fractions <- function(cytof_wide, lab_ids, reference) {
  stopifnot(reference %in% colnames(cytof_wide), all(lab_ids %in% rownames(cytof_wide)))
  cytof_wide[lab_ids, setdiff(colnames(cytof_wide), reference), drop = FALSE]
}

#' Syntactic, unique column names for cell fractions. Cluster names that differ
#' only by a sign (CD27+ B vs CD27- B) would otherwise collapse to one name and
#' the second would be silently dropped from the model.
fraction_colnames <- function(fraction_names) make.names(fraction_names, unique = TRUE)

#' DESeq2 design with karyotype as the last term. Fraction names go through
#' fraction_colnames() so they can be column names of colData.
design_formula <- function(covariates, fraction_names) {
  terms <- c(covariates, fraction_colnames(fraction_names), "karyotype")
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

#' Remove covariate contributions from every row of Y while keeping karyotype.
#'
#' Fits Y ~ 1 + karyotype + X by least squares for all genes at once and
#' subtracts X %*% beta_X. Karyotype is in the fit so its effect is not
#' absorbed by covariates that differ between groups, and it is left in the
#' result. With no covariates Y is returned unchanged.
adjust_expression <- function(Y, karyotype, X) {
  if (is.null(X) || ncol(as.matrix(X)) == 0) return(Y)
  X <- as.matrix(X)
  stopifnot(is.matrix(Y), nrow(X) == ncol(Y), length(karyotype) == ncol(Y))
  D <- cbind(intercept = 1, t21 = as.numeric(karyotype), X)
  B <- coef_rows(Y, D)
  Y - B[, -(1:2), drop = FALSE] %*% t(X)
}

#' Decompose the change in the karyotype coefficient into per-covariate terms.
#'
#' For least squares, (unadjusted - adjusted) karyotype coefficient equals
#' sum_k beta_k * (mean_k(T21) - mean_k(control)), where beta_k is covariate
#' k's coefficient in the adjusted fit. Descriptive only: it says which
#' covariates carry the difference between the two estimates, not why.
#'
#' @return data.frame(term, coef, group_diff, contribution)
attribution <- function(y, karyotype, X) {
  X <- as.matrix(X); g <- as.logical(karyotype)
  stopifnot(length(y) == nrow(X), length(g) == nrow(X))
  D <- cbind(1, as.numeric(g), X)
  b <- as.vector(coef_rows(matrix(y, nrow = 1), D))[-(1:2)]
  d <- colMeans(X[g, , drop = FALSE]) - colMeans(X[!g, , drop = FALSE])
  data.frame(term = colnames(X), coef = b, group_diff = unname(d),
             contribution = b * unname(d), stringsAsFactors = FALSE)
}
