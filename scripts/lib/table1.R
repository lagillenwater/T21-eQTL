# table1.R
#
# Summary helpers for the analysis-cohort characteristics table.
#
# Continuous variables are median [IQR] with a rank-based test, not mean (SD)
# with a t-test: Age_at_visit and BMI are right-skewed. Categorical comparisons
# use Fisher exact rather than chi-square so the helpers stay valid on small
# strata.

summarize_continuous <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
  sprintf("%.*f [%.*f-%.*f]", digits, q[2], digits, q[1], digits, q[3])
}

summarize_categorical <- function(x, level) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  n <- sum(as.character(x) == as.character(level))
  sprintf("%d (%.1f%%)", n, 100 * n / length(x))
}

compare_groups <- function(x, group) {
  g <- factor(group)
  if (nlevels(g) != 2) stop("compare_groups needs exactly two groups; got ", nlevels(g))
  keep <- !is.na(x)
  x <- x[keep]; g <- droplevels(g[keep])
  if (nlevels(g) != 2) return(list(p = NA_real_, test = "not comparable"))
  if (is.numeric(x)) {
    list(p = suppressWarnings(wilcox.test(x ~ g)$p.value), test = "Wilcoxon rank-sum")
  } else {
    tab <- table(as.character(x), g)
    list(p = tryCatch(fisher.test(tab, simulate.p.value = nrow(tab) > 2)$p.value,
                      error = function(e) NA_real_),
         test = "Fisher exact")
  }
}
