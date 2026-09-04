# sankey_flow.R
#
# The chr21 lane flow (Chr21 protein-coding -> classification -> sub-category
# -> eQTL terminal) and its SankeyMATIC serialisation, extracted from
# scripts/05_sankeymatic_export.R so they can be unit-tested. The in-script
# version broke silently on 2026-09-01: level 3 was relabelled to tier
# notation while the level-4 rule still matched the old strings, so every DE
# gene fell through to NA and the export lost both DE arms and all three eQTL
# terminals. tests/testthat/test-sankey-flow.R pins the flows drawn in
# docs/figures/Sankey.png.
#
# Node names are the short labels used in that figure. The three eQTL
# terminals are shared by the DE high and DE low arms.

FLOW_ROOT <- "Chr21 protein-coding"

LEVEL2_ORDER <- c("Outside dosage expectation", "Not assessable",
                  "Expected dosage")
LEVEL3_ORDER <- c("DE high", "Not DE", "DE low",
                  "High repeats", "Low expression", "Not estimable",
                  "Expected dosage")
LEVEL4_ORDER <- c("cis eQTL", "no cis eQTL", "no GTEx QTL", "Not DE",
                  "High repeats", "Low expression", "Not estimable",
                  "Expected dosage")
DE_LEVEL3 <- c("DE high", "DE low")

#' Aggregate the lane table into a level2 / level3 / level4 flow table.
#'
#' @param lanes data.frame or data.table with sig_lane and eqtl_lane
#'   (results/tables/chr21_lane_assignments.csv)
#' @return data.table(level2, level3, level4, n_genes), character columns,
#'   one row per distinct path, ordered top-down as the diagram reads.
#'   Non-DE lanes carry their level-3 label through to level 4 (SankeyMATIC
#'   export skips those pass-throughs; the flow CSV keeps them).
#' @section Errors: an unknown sig_lane, or a DE gene whose eqtl_lane has no
#'   terminal, stops with an error rather than dropping the gene.
lane_flow_table <- function(lanes) {
  missing <- setdiff(c("sig_lane", "eqtl_lane"), names(lanes))
  if (length(missing)) {
    stop("lanes is missing column(s): ", paste(missing, collapse = ", "))
  }
  d <- data.table(sig_lane  = as.character(lanes$sig_lane),
                  eqtl_lane = as.character(lanes$eqtl_lane))

  d[, level2 := fcase(
    sig_lane == "Expected_dosage",
      "Expected dosage",
    sig_lane %in% c("High_repeats", "Low_expression", "Not_assessable"),
      "Not assessable",
    sig_lane %in% c("DE_high", "DE_low", "Not_DE_outside_noise"),
      "Outside dosage expectation",
    default = NA_character_)]
  d[, level3 := fcase(
    sig_lane == "Expected_dosage",      "Expected dosage",
    sig_lane == "High_repeats",         "High repeats",
    sig_lane == "Low_expression",       "Low expression",
    sig_lane == "Not_assessable",       "Not estimable",
    sig_lane == "DE_high",              "DE high",
    sig_lane == "DE_low",               "DE low",
    sig_lane == "Not_DE_outside_noise", "Not DE",
    default = NA_character_)]
  if (anyNA(d$level3)) {
    stop("unknown sig_lane value(s): ",
         paste(unique(d[is.na(level3), sig_lane]), collapse = ", "))
  }

  d[, level4 := fcase(
    level3 %in% DE_LEVEL3 & eqtl_lane == "cis_eqtl",     "cis eQTL",
    level3 %in% DE_LEVEL3 & eqtl_lane == "no_cis_eqtl",  "no cis eQTL",
    level3 %in% DE_LEVEL3 & eqtl_lane == "no_GTEx_data", "no GTEx QTL",
    !level3 %in% DE_LEVEL3,                              level3,
    default = NA_character_)]
  if (anyNA(d$level4)) {
    stop("DE gene(s) with no eQTL terminal; unknown eqtl_lane value(s): ",
         paste(unique(d[is.na(level4), eqtl_lane]), collapse = ", "))
  }

  flow <- d[, .(n_genes = .N), by = .(level2, level3, level4)]
  flow[, `:=`(o2 = match(level2, LEVEL2_ORDER),
              o3 = match(level3, LEVEL3_ORDER),
              o4 = match(level4, LEVEL4_ORDER))]
  setorder(flow, o2, o3, o4)
  flow[, c("o2", "o3", "o4") := NULL]
  flow[]
}

#' Serialise a flow table as SankeyMATIC "Source [n] Target" lines.
#'
#' @param flow output of lane_flow_table()
#' @param root name of the level-1 node
#' @return character vector, one edge per line. Level-to-level pass-throughs
#'   (source name == target name) are skipped: SankeyMATIC would draw them
#'   as self-loops.
sankeymatic_lines <- function(flow, root = FLOW_ROOT) {
  l12 <- flow[, .(n = sum(n_genes)), by = .(target = level2)][, source := root]
  l23 <- flow[level2 != level3, .(n = sum(n_genes)),
              by = .(source = level2, target = level3)]
  l34 <- flow[level3 != level4, .(n = sum(n_genes)),
              by = .(source = level3, target = level4)]
  edges <- rbindlist(list(l12, l23, l34), use.names = TRUE)
  sprintf("%s [%d] %s", edges$source, edges$n, edges$target)
}
