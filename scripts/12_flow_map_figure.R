# 12_flow_map_figure.R
#
# Purpose: One two-panel figure combining the lane-flow Sankey (A) with the
#          chr21 deviating-gene map (B), for the manuscript.
#
#   A  The Sankey rendered at https://sankeymatic.com/build/ from script 05's
#      text export. That render is a manual step, so this script reads the
#      tracked PNG rather than drawing it. Both panels are raster for that
#      reason; the map is read back from script 11's PNG so the two match in
#      resolution.
#   B  chr21_deviating_map, written by script 11 for this run.
#
# The Sankey must be the render of THIS run's export, which nothing can check
# automatically - the file is produced by hand on a website. The script prints
# the file it used and when it was last written, so a stale panel A is visible
# rather than silent. Prefer a run-prefixed render (docs/figures/<run>_Sankey.png);
# the unprefixed docs/figures/Sankey.png is the historical baseline one.
#
# Inputs:  docs/figures/<run>_Sankey.png (or --sankey <path>)
#          results/runs/<run>/figures/chr21_deviating_map.png
# Outputs: results/runs/<run>/figures/chr21_flow_and_map.{pdf,png}
#
# Usage:   Rscript scripts/12_flow_map_figure.R --run adjusted
#          Rscript scripts/12_flow_map_figure.R --run adjusted --sankey path/to.png

suppressPackageStartupMessages({
  library(png)
  library(grid)
})
source("scripts/lib/run.R"); run <- load_run()

args <- commandArgs(trailingOnly = TRUE)
i    <- match("--sankey", args)
sankey_path <- if (!is.na(i)) {
  if (length(args) <= i) stop("--sankey given without a file", call. = FALSE)
  args[i + 1]
} else {
  prefixed <- file.path("docs", "figures", sprintf("%s_Sankey.png", run$name))
  if (file.exists(prefixed)) prefixed else file.path("docs", "figures", "Sankey.png")
}
map_path <- run$figure("chr21_deviating_map.png")

cat(sprintf("=== T21-eQTL: lane flow + chr21 map [run: %s] ===\n\n", run$name))

for (f in c(sankey_path, map_path))
  if (!file.exists(f))
    stop("missing panel input: ", f,
         if (identical(f, map_path)) " - run scripts/11_eqtl_figures.R first" else
           " - render script 05's export at https://sankeymatic.com/build/ and save it there",
         call. = FALSE)

# The Sankey is a hand-made render, so say out loud which file went in.
cat(sprintf("  A (Sankey): %s  [modified %s]\n", sankey_path,
            format(file.mtime(sankey_path), "%Y-%m-%d %H:%M")))
cat(sprintf("  B (map):    %s  [modified %s]\n", map_path,
            format(file.mtime(map_path), "%Y-%m-%d %H:%M")))
if (!startsWith(basename(sankey_path), paste0(run$name, "_")))
  cat(sprintf("  NOTE: panel A is not a %s-prefixed render. Confirm it was made from\n",
              run$name),
      sprintf("        %s\n", run$table("chr21_lane_sankeymatic_input.txt")))

# =============================================================================
# Compose
# =============================================================================
# Panel heights follow each image's own aspect ratio at a shared width, so
# neither panel is stretched, and the output is sized to the wider image's
# pixel width at 200 dpi so nothing is upscaled.

panel_of <- function(path) {
  img <- readPNG(path)
  list(grob = rasterGrob(img, interpolate = TRUE),
       w = dim(img)[2], h = dim(img)[1])
}
a <- panel_of(sankey_path)
b <- panel_of(map_path)
cat(sprintf("  panel pixels: A %dx%d, B %dx%d\n", a$w, a$h, b$w, b$h))

DPI    <- 200
TAG_W  <- 0.3                               # in, the gutter the A/B tags live in
panel_w <- max(a$w, b$w) / DPI              # panels at their own pixel width
aspect  <- c(a$h / a$w, b$h / b$w)          # height per unit width, per panel
fig_w   <- panel_w + TAG_W
fig_h   <- panel_w * sum(aspect)

# Drawn with grid rather than patchwork: both panels carry their own title
# text hard against the left edge of the image, and a patchwork tag is placed
# inside the plot area, so it lands on top of that text. Here the tags get a
# column of their own and cannot overlap anything.
draw_fig <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    2, 2, widths = unit.c(unit(TAG_W, "in"), unit(1, "null")),
    heights = unit(aspect, "null"))))
  draw_panel <- function(row, grob, tag) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = 2))
    grid.draw(grob)
    popViewport()
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = 1))
    grid.text(tag, x = 0.55, y = unit(1, "npc") - unit(4, "pt"),
              just = c("center", "top"), gp = gpar(fontface = "bold", cex = 1.3))
    popViewport()
  }
  draw_panel(1, a$grob, "A")
  draw_panel(2, b$grob, "B")
  popViewport()
}

stem <- run$figure("chr21_flow_and_map")
png(paste0(stem, ".png"), width = fig_w, height = fig_h, units = "in", res = DPI)
draw_fig(); invisible(dev.off())
pdf(paste0(stem, ".pdf"), width = fig_w, height = fig_h)
draw_fig(); invisible(dev.off())
for (f in paste0(stem, c(".pdf", ".png"))) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat(sprintf("  Saved: %s\n", f))
}
cat(sprintf("  figure: %.1f x %.1f in\n", fig_w, fig_h))

cat("\n=== Combined figure complete ===\n")
