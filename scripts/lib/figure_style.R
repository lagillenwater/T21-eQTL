# figure_style.R
#
# Shared figure conventions for the manuscript figures. Plots carry no
# in-plot titles or subtitles, only panel letters: the manuscript captions
# say what they said. strip_titles() takes them off a panel and records the
# text (with the numbers the script computed) in a tab-separated file,
# figures/removed_titles.txt, one row per (figure, panel), so the captions
# can be checked against it. Rerunning a script replaces its rows.

#' Remove a panel's title and subtitle, recording them.
#'
#' @param p a ggplot
#' @param figure the figure's output stem (e.g. "volcano_chr21")
#' @param panel the panel letter, "" for a one-panel figure
#' @param path the record file
#' @param tag set the panel letter as the plot tag (for panels whose letter
#'   was part of the title)
#' @return p without title and subtitle
strip_titles <- function(p, figure, panel, path, tag = FALSE) {
  record_titles(figure, panel, p$labels$title, p$labels$subtitle, path)
  p <- p + ggplot2::labs(title = NULL, subtitle = NULL)
  if (tag) p <- p + ggplot2::labs(tag = panel)
  p
}

#' Record a removed title and subtitle (also for figure-level text, e.g. a
#' patchwork plot_annotation, that is not on a ggplot); panel "" for the
#' whole figure. Replaces any earlier row for the same figure and panel.
record_titles <- function(figure, panel, title, subtitle, path) {
  flat <- function(x) if (is.null(x) || identical(x, "")) "" else gsub("\n", " / ", paste(deparse_label(x), collapse = " "))
  row <- data.table::data.table(figure = figure, panel = panel, title = flat(title), subtitle = flat(subtitle))
  old <- if (file.exists(path)) data.table::fread(path, sep = "\t", colClasses = "character") else row[0]
  old <- old[!(figure == row$figure & panel == row$panel)]
  rec <- rbind(old, row)
  data.table::setorder(rec, figure, panel)
  data.table::fwrite(rec, path, sep = "\t", quote = FALSE)
  invisible(rec)
}

#' A title or subtitle as text: plain strings as they are, plotmath
#' expressions deparsed.
deparse_label <- function(x) if (is.character(x)) x else deparse(x)
