# Tests for scripts/lib/figure_style.R: titles and subtitles come off the
# plot and into a tab-separated record, one row per (figure, panel).

skip_if_not_installed("ggplot2")   # CI installs data.table and testthat only
library(data.table)
suppressWarnings(library(ggplot2))   # "built under R x.y" only

base_plot <- function(title = "A title", subtitle = "a subtitle")
  ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point() + labs(title = title, subtitle = subtitle)

test_that("strip_titles removes the title and subtitle and keeps the rest", {
  path <- tempfile(fileext = ".txt")
  p <- strip_titles(base_plot() + labs(x = "kept"), "fig1", "A", path)
  expect_null(p$labels$title)
  expect_null(p$labels$subtitle)
  expect_equal(p$labels$x, "kept")
})

test_that("strip_titles records the removed text by figure and panel", {
  path <- tempfile(fileext = ".txt")
  strip_titles(base_plot("Title A", "Sub A"), "fig1", "A", path)
  strip_titles(base_plot("Title B", NULL), "fig1", "B", path)
  got <- fread(path)
  expect_equal(got$figure, c("fig1", "fig1"))
  expect_equal(got$panel, c("A", "B"))
  expect_equal(got$title, c("Title A", "Title B"))
  expect_equal(got$subtitle, c("Sub A", ""))
})

test_that("strip_titles replaces a figure's earlier row instead of adding a second", {
  path <- tempfile(fileext = ".txt")
  strip_titles(base_plot("old"), "fig1", "A", path)
  strip_titles(base_plot("other figure"), "fig2", "A", path)
  strip_titles(base_plot("new"), "fig1", "A", path)
  got <- fread(path)
  expect_equal(nrow(got), 2L)
  expect_equal(got[figure == "fig1", title], "new")
})

test_that("strip_titles sets the panel letter as the tag when asked", {
  p <- strip_titles(base_plot(), "fig1", "C", tempfile(), tag = TRUE)
  expect_equal(p$labels$tag, "C")
})

test_that("strip_titles keeps line breaks of a two-line subtitle readable", {
  path <- tempfile(fileext = ".txt")
  strip_titles(base_plot("T", "line one\nline two"), "fig1", "A", path)
  expect_equal(fread(path)$subtitle, "line one / line two")
})
