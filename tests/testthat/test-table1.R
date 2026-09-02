
test_that("summarize_continuous formats median and IQR", {
  expect_equal(summarize_continuous(c(1,2,3,4,5,6,7,8,9,10), digits = 1), "5.5 [3.2-7.8]")
})

test_that("summarize_continuous excludes missing values", {
  # median [IQR] of the 4 non-missing values (1,2,3,5), type-7 quantile:
  # median = 2.5, Q1 = 1.75 -> "1.8", Q3 = 3.5. Verified by hand against
  # quantile(c(1,2,3,5), c(.25,.5,.75), type = 7).
  expect_match(summarize_continuous(c(1, 2, 3, NA, 5), digits = 1), "^2\\.5 ")
})

test_that("summarize_continuous handles an all-missing vector", {
  expect_equal(summarize_continuous(c(NA_real_, NA_real_)), "-")
})

test_that("summarize_categorical gives count and percent of non-missing", {
  x <- c("Male", "Female", "Female", "Female", NA)
  expect_equal(summarize_categorical(x, "Female"), "3 (75.0%)")
  expect_equal(summarize_categorical(x, "Male"), "1 (25.0%)")
})

test_that("summarize_categorical handles an absent level", {
  expect_equal(summarize_categorical(c("A", "A"), "B"), "0 (0.0%)")
})

test_that("compare_groups uses Wilcoxon for continuous and Fisher for categorical", {
  set.seed(20)
  cont <- compare_groups(c(rnorm(30), rnorm(30, 3)), rep(c("a", "b"), each = 30))
  expect_equal(cont$test, "Wilcoxon rank-sum")
  expect_lt(cont$p, 0.001)
  cat2 <- compare_groups(c(rep("M", 25), rep("F", 5), rep("M", 5), rep("F", 25)),
                         rep(c("a", "b"), each = 30))
  expect_equal(cat2$test, "Fisher exact")
  expect_lt(cat2$p, 0.001)
})

test_that("compare_groups refuses anything but two groups", {
  expect_error(compare_groups(1:10, rep("a", 10)), "exactly two groups")
})
