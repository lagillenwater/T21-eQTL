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

#' Build the age/BMI/sex/sample-source/study-visit/comorbidity rows of
#' Table 1 for a given cohort. Shared by scripts/00_preprocess_data.R (the
#' pre-run analysis cohort) and scripts/S2_table1_run_cohort.R (a run's
#' post-exclusion cohort) so both tables are built the same way; only the
#' cohort passed in differs. Does not include cohort-size header rows -
#' callers add those, since what counts as "excluded" differs between them.
#'
#' @param cohort data.table/data.frame with Karyotype ("T21"/"Control"),
#'   RecordID, Age_at_visit, BMI, Sex, Sample_source, Event_name.
#' @param como_path optional comorbidity TSV
#'   (data/P4C_Comorbidity_020921.tsv format); the block is skipped if NULL,
#'   the file is absent, or its columns don't match.
#' @return data.table of Table 1 rows.
build_table1_rows <- function(cohort, como_path = NULL) {
  cohort <- data.table::as.data.table(cohort)
  t21a <- cohort[Karyotype == "T21"]
  ctla <- cohort[Karyotype == "Control"]

  row_cont <- function(label, var, digits = 1) {
    cmp <- compare_groups(cohort[[var]], cohort$Karyotype)
    data.table::data.table(characteristic = label, level = "median [IQR]",
               t21 = summarize_continuous(t21a[[var]], digits),
               control = summarize_continuous(ctla[[var]], digits),
               overall = summarize_continuous(cohort[[var]], digits),
               p_value = cmp$p, test = cmp$test,
               n_missing = sum(is.na(cohort[[var]])))
  }
  row_cat <- function(label, var) {
    cmp  <- compare_groups(cohort[[var]], cohort$Karyotype)
    levs <- sort(unique(as.character(cohort[[var]][!is.na(cohort[[var]])])))
    data.table::rbindlist(lapply(seq_along(levs), function(i)
      data.table::data.table(characteristic = if (i == 1) label else "", level = levs[i],
                 t21 = summarize_categorical(t21a[[var]], levs[i]),
                 control = summarize_categorical(ctla[[var]], levs[i]),
                 overall = summarize_categorical(cohort[[var]], levs[i]),
                 p_value = if (i == 1) cmp$p else NA_real_,
                 test = if (i == 1) cmp$test else "",
                 n_missing = if (i == 1) sum(is.na(cohort[[var]])) else NA_integer_)))
  }

  tbl1 <- data.table::rbindlist(list(
    row_cont("Age at visit (years)", "Age_at_visit"),
    row_cont("BMI", "BMI"),
    row_cat("Sex", "Sex"),
    row_cat("Sample source", "Sample_source"),
    row_cat("Study visit", "Event_name")))

  if (!is.null(como_path) && file.exists(como_path)) {
    como <- data.table::fread(como_path, skip = 1)
    if (all(c("RecordID", "Condition", "HasCondition") %in% names(como))) {
      como <- como[RecordID %in% cohort$RecordID]
      top  <- como[, .(n_with = sum(HasCondition == 1, na.rm = TRUE)), by = Condition][
        order(-n_with)][seq_len(min(5, .N))]
      wide <- data.table::dcast(como, RecordID ~ Condition, value.var = "HasCondition")
      ann  <- merge(cohort[, .(RecordID, Karyotype)], wide, by = "RecordID")
      tbl1 <- data.table::rbindlist(list(tbl1, data.table::rbindlist(lapply(top$Condition, function(cn) {
        cmp <- compare_groups(as.character(ann[[cn]]), ann$Karyotype)
        data.table::data.table(characteristic = cn, level = "n (%) with condition",
                   t21 = summarize_categorical(ann[Karyotype == "T21"][[cn]], 1),
                   control = summarize_categorical(ann[Karyotype == "Control"][[cn]], 1),
                   overall = summarize_categorical(ann[[cn]], 1),
                   p_value = cmp$p, test = cmp$test, n_missing = sum(is.na(ann[[cn]])))
      }))))
    }
  }
  tbl1
}
