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

#' DESeq2 design with karyotype after the main effects. Fraction names go
#' through fraction_colnames() so they can be column names of colData.
#'
#' With `interact_karyotype`, a karyotype x fraction term follows for every
#' fraction, so each cell population may relate to expression differently in
#' T21 and D21. The karyotype coefficient is then the T21 effect where every
#' fraction is 0; script 01 centres the fractions, so that is the cohort-mean
#' composition, and its contrast name is unchanged.
design_formula <- function(covariates, fraction_names, interact_karyotype = FALSE) {
  fr <- fraction_colnames(fraction_names)
  terms <- c(covariates, fr, "karyotype")
  if (interact_karyotype && length(fr)) terms <- c(terms, paste0("karyotype:", fr))
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

#' Remove covariate contributions from every row of Y while keeping karyotype.
#'
#' Fits Y ~ 1 + karyotype + X by least squares for all genes at once and
#' subtracts X %*% beta_X. Karyotype is in the fit so its effect is not
#' absorbed by covariates that differ between groups, and it is left in the
#' result. With no covariates Y is returned unchanged.
#'
#' `interact` names columns of X whose slope may differ by karyotype: their
#' karyotype x covariate products join the fit and are subtracted too, so no
#' covariate effect remains within either group. The karyotype effect kept is
#' the one where those columns are 0 (centre them for the mean composition).
adjust_expression <- function(Y, karyotype, X, interact = NULL) {
  if (is.null(X) || ncol(as.matrix(X)) == 0) return(Y)
  X <- as.matrix(X)
  stopifnot(is.matrix(Y), nrow(X) == ncol(Y), length(karyotype) == ncol(Y))
  missing_cols <- setdiff(interact, colnames(X))
  if (length(missing_cols))
    stop("interaction column(s) not among the covariates: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  t21 <- as.numeric(karyotype)
  if (length(interact)) {
    XI <- X[, interact, drop = FALSE] * t21
    colnames(XI) <- paste0("t21:", interact)
    X <- cbind(X, XI)
  }
  D <- cbind(intercept = 1, t21 = t21, X)
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

#' The run's covariates as a numeric matrix (no intercept), one row per row of
#' `cohort`: model.matrix() dummies for the metadata covariates
#' (run$covariates) and, in composition runs, the cell fractions of
#' processed/composition_fractions.csv matched by LabID. The same covariates
#' script 01 removes from expression; the within-T21 fits project them out of
#' genotype too. NULL when the run has none.
run_covariate_matrix <- function(run, cohort) {
  cohort <- as.data.frame(cohort)
  X <- NULL
  if (length(run$covariates)) {
    vals <- cohort[, run$covariates, drop = FALSE]
    if (anyNA(vals)) stop("missing covariate values in the cohort", call. = FALSE)
    X <- stats::model.matrix(~ ., vals)[, -1, drop = FALSE]
  }
  if (!is.null(run$composition)) {
    comp <- data.table::fread(run$processed("composition_fractions.csv"))
    rows <- match(cohort$LabID, comp$LabID)
    if (anyNA(rows)) stop("cohort samples without composition fractions", call. = FALSE)
    X <- cbind(X, as.matrix(comp[rows, setdiff(names(comp), "LabID"), with = FALSE]))
  }
  if (!is.null(X) && anyNA(X)) stop("missing covariate values in the cohort", call. = FALSE)
  X
}

#' Display labels for run covariates (columns of run_covariate_matrix): plain
#' names for the metadata dummies, and Greek letters spelled out, because the
#' default pdf() device cannot draw them and cairo is not available on every
#' machine here (scripts/07_three_panel_figure.R).
display_covariate <- function(x) {
  x <- sub("^Age_at_visit$", "Age", x)
  x <- sub("^SexMale$", "Sex (male)", x)
  x <- sub("^Sample_source", "Source: ", x)
  gsub("\u03b3\u03b4", "gamma-delta", x)
}
