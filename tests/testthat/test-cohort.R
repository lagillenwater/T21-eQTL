
test_that("subject_id_from_labid strips the trailing visit suffix", {
  expect_equal(subject_id_from_labid(c("HTP0001A", "HTP0002B2", "HTP0003")),
               c("HTP0001", "HTP0002", "HTP0003"))
})

test_that("wgs_subjects reads headers only and drops the fixed VCF columns", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,HTP0001A,HTP0002B",
               "chr21,100,.,A,G,.,PASS,.,GT,0/1,1/1"), f)
  expect_setequal(wgs_subjects(f), c("HTP0001", "HTP0002"))
})

test_that("wgs_subjects unions across files and de-duplicates", {
  f1 <- tempfile(fileext = ".csv"); f2 <- tempfile(fileext = ".csv")
  writeLines("CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,HTP0001A", f1)
  writeLines("CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,HTP0001B,HTP0009A", f2)
  expect_setequal(wgs_subjects(c(f1, f2)), c("HTP0001", "HTP0009"))
})

test_that("analysis_cohort keeps genotyped T21 and ALL controls", {
  meta <- data.table::data.table(
    LabID     = c("A1A", "A2A", "A3A", "C1A", "C2A"),
    Karyotype = c("T21", "T21", "T21", "Control", "Control"),
    has_wgs   = c(TRUE, TRUE, FALSE, TRUE, FALSE))
  got <- analysis_cohort(meta)
  expect_equal(sum(got$Karyotype == "T21"), 2L)      # A3A dropped, no WGS
  expect_equal(sum(got$Karyotype == "Control"), 2L)  # both kept regardless
})

test_that("analysis_cohort refuses metadata missing has_wgs", {
  meta <- data.table::data.table(LabID = "A1A", Karyotype = "T21")
  expect_error(analysis_cohort(meta), "has_wgs")
})

test_that("build_cohort applies keep rule, exclusions, covariate completeness and cytof match with reasons", {
  meta <- data.table::data.table(
    LabID = c("T1A", "T2A", "T3A", "T4A", "C1A", "C2A"),
    Karyotype = c("T21", "T21", "T21", "T21", "Control", "Control"),
    has_wgs = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
    karyotype_subtype = c("T21", "mosaic_T21", "T21", "T21", "D21", "D21"),
    Age_at_visit = c(20, 21, 22, 23, 24, 25), BMI = c(25, 26, 27, NA, 28, 29))
  run <- list(covariates = c("Age_at_visit", "BMI"),
              composition = list(source = "cytof"),
              exclude = list(karyotype_subtype = "mosaic_T21"))
  got <- build_cohort(meta, run, cytof_ids = c("T1A", "T2A", "T3A", "T4A", "C1A"))
  expect_equal(nrow(got), 6)
  expect_equal(got$excluded_reason,
               c(NA, "exclude:karyotype_subtype=mosaic_T21", "keep_rule:no_wgs",
                 "missing:BMI", NA, "missing:cytof"))
  base <- build_cohort(meta, list(covariates = character(0), composition = NULL, exclude = list()))
  expect_equal(is.na(base$excluded_reason), c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
})
