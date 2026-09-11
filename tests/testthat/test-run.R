base_cfg <- function() list(
  name = "t", covariates = character(0), composition = NULL, exclude = list(),
  ploidy = 1.5,
  thresholds = list(alpha_de = 0.01, deviation_lfc = log2(1.5),
                    deviation_lfc_t2 = log2(4/3), low_expr_basemean = 30,
                    fdr_gene = 0.05, gtex_pval_keep = 1e-4,
                    outlier_fdr = 0.10, alpha_repro = 0.05, n_positive_controls = 10, decoy_min_distance = 5e6))

cfg_text <- function(name) sprintf(
  'list(name = "%s", covariates = character(0), composition = NULL, exclude = list(), ploidy = 1.5, thresholds = list(alpha_de = 0.01, deviation_lfc = log2(1.5), deviation_lfc_t2 = log2(4/3), low_expr_basemean = 30, fdr_gene = 0.05, gtex_pval_keep = 1e-4, outlier_fdr = 0.10, alpha_repro = 0.05, n_positive_controls = 10, decoy_min_distance = 5e6))',
  name)

make_root <- function(names) {
  root <- tempfile(); dir.create(file.path(root, "config", "runs"), recursive = TRUE)
  for (n in names) writeLines(cfg_text(n), file.path(root, "config", "runs", paste0(n, ".R")))
  root
}

test_that("validate_run_config accepts the baseline shape and rejects deviations", {
  expect_identical(validate_run_config(base_cfg()), base_cfg())
  bad <- base_cfg(); bad$thresholds$fdr_gene <- NULL
  expect_error(validate_run_config(bad), "run config: missing threshold fdr_gene")
  bad <- base_cfg(); bad$covariates <- c("Age_at_visit", "Height")
  expect_error(validate_run_config(bad), "run config: unknown covariate Height")
  bad <- base_cfg(); bad$exclude <- list(mosaic_T21 = TRUE)
  expect_error(validate_run_config(bad), "run config: exclude values must be character")
  bad <- base_cfg(); bad$name <- "has space"
  expect_error(validate_run_config(bad), "run config: name must match")
  bad <- base_cfg(); bad$thresholds$extra <- 1
  expect_error(validate_run_config(bad), "run config: unknown threshold extra")
})

test_that("load_run parses --run, defaults to baseline, and creates run dirs", {
  root <- make_root("demo")
  run <- load_run(c("--run", "demo"), root = root)
  expect_equal(run$name, "demo")
  expect_true(dir.exists(file.path(root, "results", "runs", "demo", "tables")))
  expect_equal(run$table("x.csv"), file.path(root, "results", "runs", "demo", "tables", "x.csv"))
  expect_equal(run$figure("f"), file.path(root, "results", "runs", "demo", "figures", "f"))
  expect_equal(run$processed("p.csv"), file.path(root, "results", "runs", "demo", "processed", "p.csv"))
  expect_error(load_run(c("--run", "nope"), root = root), "run config: no file")
  expect_error(load_run(character(0), root = root), "run config: no file.*baseline")
  expect_error(load_run(c("--run"), root = root), "without a run name")
  expect_error(load_run(c("--run", ""), root = root), "without a run name")
})

test_that("T21_RUN is honoured when --run is absent, and --run wins over it", {
  root <- make_root(c("a", "b"))
  old <- Sys.getenv("T21_RUN", unset = NA)
  Sys.setenv(T21_RUN = "a")
  on.exit(if (is.na(old)) Sys.unsetenv("T21_RUN") else Sys.setenv(T21_RUN = old))
  expect_equal(load_run(character(0), root = root)$name, "a")
  expect_equal(load_run(c("--run", "b"), root = root)$name, "b")
})

test_that("run$partner_pool applies the mean-count and non-chr21 rules", {
  root <- make_root("demo")
  run <- load_run(c("--run", "demo"), root = root)
  counts <- data.table::data.table(EnsemblID = c("e1", "e2", "e3"), Gene_name = c("A", "B", "C"),
                                   Chr = c("chr1", "chr21", "chr2"),
                                   s1 = c(100, 100, 1), s2 = c(100, 100, 1))
  expect_equal(run$partner_pool(counts), c(TRUE, FALSE, FALSE))
  # restricting to samples changes the mean: gene C has 1 in s1, 100 in s3
  counts$s3 <- c(100, 100, 100)
  expect_equal(run$partner_pool(counts, lab_ids = c("s1", "s2")), c(TRUE, FALSE, FALSE))
  expect_equal(run$partner_pool(counts, lab_ids = c("s1", "s3")), c(TRUE, FALSE, TRUE))
})

test_that("run$expression reads the artifact as a genes x samples matrix", {
  root <- make_root("demo")
  run <- load_run(c("--run", "demo"), root = root)
  data.table::fwrite(data.table::data.table(Gene_name = c("A", "B"), s1 = c(1, 2), s2 = c(3, 4)),
                     run$processed("expression_adjusted.csv"))
  E <- run$expression()
  expect_equal(dim(E), c(2, 2))
  expect_equal(rownames(E), c("A", "B"))
  expect_equal(unname(E["B", "s2"]), 4)
})
