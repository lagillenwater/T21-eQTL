# 13_alluvial_lane_assignment.R
#
# Purpose: Comprehensive Sankey/alluvial diagram for chr21 genes that mirrors
#          the paper's Panel D (Hunter et al. 2023, BMC Biology 21:228) but
#          extends it symmetrically to include both arms (DE_low and DE_high)
#          and adds a "no_GTEx_data" terminal so genes that cannot be
#          evaluated for an eQTL explanation are not silently miscategorized.
#
#          Flow structure:
#            Level 1: All chr21 genes
#            Level 2: raw_FC arm (>= 1.5 vs < 1.5)
#            Level 3: categorization (Expected dosage / High repeats /
#                     Low expression / Not DE / DE)
#            Level 4: eQTL terminal for DE arms
#                     (Explained / Unexplained / No GTEx data); other
#                     categories pass through unchanged.
#
# Inputs:
#   - results/tables/chr21_lane_assignments.csv
#
# Outputs:
#   - results/figures/chr21_lane_alluvial.pdf
#   - results/figures/chr21_lane_alluvial.png
#   - results/tables/chr21_lane_alluvial_flow.csv
#   - results/figures/chr21_lane_alluvial_session_info.txt
#
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
})

set.seed(42)

cat("=== T21-eQTL: Comprehensive Lane Alluvial ===\n\n")

# =============================================================================
# STEP 1: Load lane assignments
# =============================================================================

cat("Step 1: Loading lane assignments...\n")

lanes <- fread("results/tables/chr21_lane_assignments.csv")
cat(sprintf("  chr21 genes: %d\n", nrow(lanes)))

# =============================================================================
# STEP 2: Build flow data - cohort-noise filter is the FIRST split
# =============================================================================

cat("\nStep 2: Building flow structure...\n")

# Level 2: cohort-noise filter (first categorization). Genes that fail are
# categorized as Expected dosage; genes that pass continue to further
# categorization downstream.
lanes[, level2 := fifelse(passes_magnitude_filter,
                          "Outside cohort noise",
                          "Within cohort noise (Expected dosage)")]

# Level 3: only meaningful for "Outside cohort noise" survivors.
#   For survivors: split into DE_high / DE_low / High_repeats /
#                  Low_expression / Not_DE_outside_noise.
#   For Expected dosage: terminate at level 2 (pass through unchanged).
lanes[, level3 := fcase(
  level2 == "Within cohort noise (Expected dosage)",
    "Within cohort noise (Expected dosage)",
  sig_lane == "DE_low",
    "DE (low) [< 1.5 raw FC]",
  sig_lane == "DE_high",
    "DE (high) [>= 1.5 raw FC]",
  sig_lane == "High_repeats",
    "High repeats",
  sig_lane == "Low_expression",
    "Low expression",
  sig_lane == "Not_DE_outside_noise",
    "Not DE (outside cohort noise)",
  default = "Other")]

# Level 4: eQTL terminal (only for DE_low / DE_high; other lanes pass through)
lanes[, level4 := fcase(
  level3 == "DE (low) [< 1.5 raw FC]" & eqtl_lane == "explained",
    "DE (low): eQTL-supported",
  level3 == "DE (low) [< 1.5 raw FC]" & eqtl_lane == "unexplained",
    "DE (low): eQTL-tested, not supported",
  level3 == "DE (low) [< 1.5 raw FC]" & eqtl_lane == "no_GTEx_data",
    "DE (low): no GTEx eQTL data",
  level3 == "DE (high) [>= 1.5 raw FC]" & eqtl_lane == "explained",
    "DE (high): eQTL-supported",
  level3 == "DE (high) [>= 1.5 raw FC]" & eqtl_lane == "unexplained",
    "DE (high): eQTL-tested, not supported",
  level3 == "DE (high) [>= 1.5 raw FC]" & eqtl_lane == "no_GTEx_data",
    "DE (high): no GTEx eQTL data",
  default = level3)]

# Aggregate counts
flow <- lanes[, .(n_genes = .N),
              by = .(level2, level3, level4)]

# Order strata for vertical layout
level2_order <- c("Outside cohort noise",
                  "Within cohort noise (Expected dosage)")
level3_order <- c("DE (high) [>= 1.5 raw FC]",
                  "High repeats", "Low expression",
                  "Not DE (outside cohort noise)",
                  "DE (low) [< 1.5 raw FC]",
                  "Within cohort noise (Expected dosage)")
level4_order <- c(
  "DE (high): eQTL-supported",
  "DE (high): eQTL-tested, not supported",
  "DE (high): no GTEx eQTL data",
  "High repeats", "Low expression",
  "Not DE (outside cohort noise)",
  "DE (low): eQTL-supported",
  "DE (low): eQTL-tested, not supported",
  "DE (low): no GTEx eQTL data",
  "Within cohort noise (Expected dosage)"
)

flow[, level2 := factor(level2, levels = level2_order)]
flow[, level3 := factor(level3, levels = level3_order)]
flow[, level4 := factor(level4, levels = level4_order)]
setorder(flow, level2, level3, level4)

fwrite(flow, "results/tables/chr21_lane_alluvial_flow.csv")
cat(sprintf("  Wrote flow table (%d rows)\n", nrow(flow)))
print(flow)

# =============================================================================
# STEP 3: Render alluvial
# =============================================================================

cat("\nStep 3: Rendering alluvial...\n")

# Long format with row id (alluvium) and a stable "final destination"
flow[, alluvium := .I]
to_long <- melt(flow, id.vars = c("alluvium", "n_genes", "level4"),
                measure.vars = c("level2", "level3", "level4"),
                variable.name = "x_axis", value.name = "stratum")
to_long[, x_axis := factor(x_axis, levels = c("level2", "level3", "level4"))]
to_long[, final_category := level4]

# Color palette - flows colored by FINAL terminal so the eye traces
# explanatory vs unexplained vs no-data outcomes through the diagram.
palette_terminal <- c(
  # Cohort-noise filter strata (level 2)
  "Outside cohort noise"                       = "#808080",
  "Within cohort noise (Expected dosage)"      = "#5AB4AC",
  # Sub-categories (level 3)
  "DE (high) [>= 1.5 raw FC]"                  = "#1F77B4",
  "DE (low) [< 1.5 raw FC]"                    = "#4DAF4A",
  "High repeats"                               = "#E8967A",
  "Low expression"                             = "#DDA0DD",
  "Not DE (outside cohort noise)"              = "#D2B48C",
  # eQTL terminals (level 4)
  "DE (high): eQTL-supported"                  = "#1F77B4",
  "DE (high): eQTL-tested, not supported"      = "#9467BD",
  "DE (high): no GTEx eQTL data"               = "#AEC7E8",
  "DE (low): eQTL-supported"                   = "#4DAF4A",
  "DE (low): eQTL-tested, not supported"       = "#FF7F0E",
  "DE (low): no GTEx eQTL data"                = "#999999"
)

p <- ggplot(to_long,
            aes(x = x_axis, stratum = stratum, alluvium = alluvium,
                y = n_genes)) +
  geom_flow(aes(fill = final_category),
            alpha = 0.7, curve_type = "cubic", segments = 100,
            width = 1/20) +
  geom_stratum(aes(fill = stratum),
               width = 1/20, color = "black", linewidth = 0.4) +
  geom_text(stat = "stratum",
            aes(label = sprintf("%s\n(%d)", after_stat(stratum),
                                after_stat(count))),
            size = 2.4, lineheight = 0.85) +
  scale_fill_manual(values = palette_terminal) +
  scale_x_discrete(limits = c("level2", "level3", "level4"),
                   labels = c("Cohort-noise filter",
                              "Sub-category", "eQTL terminal"),
                   expand = c(0.08, 0.08)) +
  labs(
    title    = "Chr21 gene classification by lane",
    subtitle = sprintf(
      paste0("Cohort-noise filter (1.0 SD of non-chr21 distribution) is ",
             "the first split: small-magnitude deviations -> Expected ",
             "dosage. Survivors then split by raw FC sign and eQTL ",
             "evidence. (n=%d)"),
      nrow(lanes)),
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9),
    axis.text.x      = element_text(size = 10, face = "bold"),
    axis.title.x     = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "none"
  )

ggsave("results/figures/chr21_lane_alluvial.pdf", p,
       width = 13, height = 9, bg = "white")
cat("  Saved: results/figures/chr21_lane_alluvial.pdf\n")

ggsave("results/figures/chr21_lane_alluvial.png", p,
       width = 13, height = 9, dpi = 150, bg = "white")
cat("  Saved: results/figures/chr21_lane_alluvial.png\n")

# =============================================================================
# STEP 4: SankeyMATIC export
# =============================================================================

cat("\nStep 4: Writing SankeyMATIC input...\n")

# Aggregate flows at each level transition into Source [count] Target lines
sm_lines <- character(0)

# Level 1 -> Level 2: All chr21 -> raw FC arms
l12 <- flow[, .(n_genes = sum(n_genes)), by = level2]
for (i in seq_len(nrow(l12))) {
  sm_lines <- c(sm_lines,
                sprintf("Chr21 protein-coding [%d] %s",
                        l12$n_genes[i], l12$level2[i]))
}

# Level 2 -> Level 3 (skip pass-throughs where level3 == level2; those
# would render as self-loops in SankeyMATIC)
l23 <- flow[as.character(level2) != as.character(level3),
            .(n_genes = sum(n_genes)),
            by = .(level2, level3)]
for (i in seq_len(nrow(l23))) {
  sm_lines <- c(sm_lines,
                sprintf("%s [%d] %s",
                        l23$level2[i], l23$n_genes[i], l23$level3[i]))
}

# Level 3 -> Level 4 (skip pass-throughs)
l34 <- flow[as.character(level3) != as.character(level4),
            .(n_genes = sum(n_genes)),
            by = .(level3, level4)]
for (i in seq_len(nrow(l34))) {
  sm_lines <- c(sm_lines,
                sprintf("%s [%d] %s",
                        l34$level3[i], l34$n_genes[i], l34$level4[i]))
}

# Color hints (as SankeyMATIC :Node #color directives in comments so they
# can be pasted into the SankeyMATIC color editor if desired). Mirrors
# palette_terminal above.
color_lines <- c(
  "",
  "// Suggested colors (paste into SankeyMATIC color editor)",
  "// :Chr21 protein-coding #808080",
  "// :Outside cohort noise #808080",
  "// :Within cohort noise (Expected dosage) #5AB4AC",
  "// :DE (high) [>= 1.5 raw FC] #1F77B4",
  "// :DE (low) [< 1.5 raw FC] #4DAF4A",
  "// :High repeats #E8967A",
  "// :Low expression #DDA0DD",
  "// :Not DE (outside cohort noise) #D2B48C",
  "// :DE (high): eQTL-supported #1F77B4",
  "// :DE (high): eQTL-tested, not supported #9467BD",
  "// :DE (high): no GTEx eQTL data #AEC7E8",
  "// :DE (low): eQTL-supported #4DAF4A",
  "// :DE (low): eQTL-tested, not supported #FF7F0E",
  "// :DE (low): no GTEx eQTL data #999999"
)

header <- c(
  paste0("// Chr21 lane assignment Sankey input"),
  paste0("// Generated by 13_alluvial_lane_assignment.R"),
  paste0("// padj < 0.01 after ploidy normalization; eQTL-supported = ",
         ">=1 cis variant matches deviation direction AND reproduces in ",
         "T21 at p < 0.05; magnitude filter at 1.0 cohort-noise SD."),
  paste0("// Total genes: ", sum(flow$n_genes)),
  ""
)

sm_path <- "results/tables/chr21_lane_sankeymatic_input.txt"
writeLines(c(header, sm_lines, color_lines), sm_path)
cat(sprintf("  Wrote %s (%d flow lines)\n", sm_path, length(sm_lines)))

# =============================================================================
# STEP 5: Verification
# =============================================================================

cat("\n=== Verification ===\n")
cat(sprintf("Total chr21 genes in flow: %d (expect %d)\n",
            sum(flow$n_genes), nrow(lanes)))
stopifnot(sum(flow$n_genes) == nrow(lanes))

cat("\nTerminal counts:\n")
print(flow[, .(n_genes = sum(n_genes)), by = level4][order(level4)])

writeLines(capture.output(sessionInfo()),
           "results/figures/chr21_lane_alluvial_session_info.txt")

cat("\n=== Alluvial complete ===\n")
