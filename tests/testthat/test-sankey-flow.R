library(data.table)

# One row per (sig_lane, eqtl_lane) combination the lane rule can produce.
lanes_all <- data.table(
  sig_lane  = c("Expected_dosage", "High_repeats", "Low_expression",
                "Not_assessable",
                "DE_high", "DE_high", "DE_high",
                "DE_low", "DE_low", "DE_low",
                "Not_DE_outside_noise"),
  eqtl_lane = c(rep("not_evaluated", 4),
                "cis_eqtl", "no_cis_eqtl", "no_GTEx_data",
                "cis_eqtl", "no_cis_eqtl", "no_GTEx_data",
                "not_evaluated")
)

# The lane composition behind docs/figures/Sankey.png (2026-09-01 run):
# 104 expected dosage, 12 high-repeat, 29 low-expression, DE high 7
# (5 cis / 1 no cis / 1 no GTEx), DE low 8 (4 / 3 / 1).
lanes_fig <- data.table(
  sig_lane  = c(rep("Expected_dosage", 104), rep("High_repeats", 12),
                rep("Low_expression", 29), rep("DE_high", 7),
                rep("DE_low", 8)),
  eqtl_lane = c(rep("not_evaluated", 145),
                rep("cis_eqtl", 5), "no_cis_eqtl", "no_GTEx_data",
                rep("cis_eqtl", 4), rep("no_cis_eqtl", 3), "no_GTEx_data")
)

parse_sankey_lines <- function(lines) {
  data.table(
    source = sub(" \\[\\d+\\] .*$", "", lines),
    n      = as.integer(sub("^.* \\[(\\d+)\\] .*$", "\\1", lines)),
    target = sub("^.* \\[\\d+\\] ", "", lines)
  )
}

test_that("every DE gene reaches an eQTL terminal", {
  # The in-script version of this mapping relabelled level 3 and left the
  # level-4 rule matching the old strings, so every DE gene fell through to
  # NA and the export lost both DE arms and all three eQTL terminals.
  flow <- lane_flow_table(lanes_all)
  de <- flow[level3 %in% c("DE high", "DE low")]
  expect_equal(nrow(de), 6)
  expect_true(all(de$level4 %in% c("cis eQTL", "no cis eQTL", "no GTEx QTL")))
  expect_false(anyNA(flow))
})

test_that("DE arms hang off Outside dosage expectation; unassessable lanes off Not assessable", {
  flow <- lane_flow_table(lanes_all)
  expect_equal(unique(flow[level3 %in% c("DE high", "DE low", "Not DE"), level2]),
               "Outside dosage expectation")
  expect_setequal(flow[level2 == "Not assessable", level3],
                  c("High repeats", "Low expression", "Not estimable"))
  expect_equal(flow[level2 == "Expected dosage", level3], "Expected dosage")
})

test_that("non-DE lanes pass through unchanged to level 4", {
  flow <- lane_flow_table(lanes_all)
  pass <- flow[!level3 %in% c("DE high", "DE low")]
  expect_equal(pass$level4, pass$level3)
})

test_that("gene count is conserved in the flow table and at every SankeyMATIC stage", {
  flow <- lane_flow_table(lanes_fig)
  expect_equal(sum(flow$n_genes), nrow(lanes_fig))
  parsed <- parse_sankey_lines(sankeymatic_lines(flow))
  expect_equal(sum(parsed[source == "Chr21 protein-coding", n]), 160)
  expect_equal(sum(parsed[source == "Outside dosage expectation", n]), 15)
  expect_equal(sum(parsed[source == "Not assessable", n]), 41)
  expect_equal(sum(parsed[source %in% c("DE high", "DE low"), n]), 15)
})

test_that("the export reproduces the flows drawn in docs/figures/Sankey.png", {
  lines <- sankeymatic_lines(lane_flow_table(lanes_fig))
  expected <- c(
    "Chr21 protein-coding [15] Outside dosage expectation",
    "Chr21 protein-coding [41] Not assessable",
    "Chr21 protein-coding [104] Expected dosage",
    "Outside dosage expectation [7] DE high",
    "Outside dosage expectation [8] DE low",
    "Not assessable [12] High repeats",
    "Not assessable [29] Low expression",
    "DE high [5] cis eQTL",
    "DE high [1] no cis eQTL",
    "DE high [1] no GTEx QTL",
    "DE low [4] cis eQTL",
    "DE low [3] no cis eQTL",
    "DE low [1] no GTEx QTL"
  )
  expect_setequal(lines, expected)
})

test_that("no line is a self-loop and the eQTL terminals are shared by both DE arms", {
  parsed <- parse_sankey_lines(sankeymatic_lines(lane_flow_table(lanes_fig)))
  expect_false(any(parsed$source == parsed$target))
  expect_equal(parsed[target == "cis eQTL", sort(source)], c("DE high", "DE low"))
  expect_equal(parsed[target == "no GTEx QTL", sort(source)], c("DE high", "DE low"))
})

test_that("an unknown sig_lane, or a DE gene without an eQTL terminal, is an error", {
  expect_error(lane_flow_table(data.table(sig_lane = "Sig_high_FC",
                                          eqtl_lane = "cis_eqtl")),
               "sig_lane")
  expect_error(lane_flow_table(data.table(sig_lane = "DE_low",
                                          eqtl_lane = "not_evaluated")),
               "eqtl_lane")
})

test_that("required columns are checked", {
  expect_error(lane_flow_table(data.table(sig_lane = "DE_low")), "eqtl_lane")
})
