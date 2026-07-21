# 06_alluvial_with_eqtl.R
#
# Purpose: Create enhanced alluvial diagram including eQTL classification
#          Extends Panel D to show which DE genes have eQTL support
#          Splits DE genes by effect size (FC < 0.75 vs FC >= 0.75)
#
# Inputs:
#   - results/tables/chr21_genes_categorized.csv
#   - results/tables/final_gene_classification.csv (with eQTL info)
#
# Outputs:
#   - results/figures/panel_D_with_eqtl.pdf
#   - results/figures/eqtl_breakdown_bar.pdf
#   - results/tables/panel_D_eqtl_summary.csv
#
# Date: 2025-12-10

# Load required libraries
library(tidyverse)
library(ggalluvial)

# Set seed for reproducibility
set.seed(42)

cat("=== T21-eQTL Analysis: Alluvial Plot with eQTL Classification ===\n\n")

# =============================================================================
# STEP 1: Load data
# =============================================================================

cat("Step 1: Loading categorized chr21 genes and eQTL results...\n")

if (!file.exists("results/tables/chr21_genes_categorized.csv")) {
  stop("Categorized data not found. Run 02_categorize_genes.R first.")
}

if (!file.exists("results/tables/final_gene_classification.csv")) {
  stop("eQTL classification not found. Run 05_eqtl_analysis.R first.")
}

chr21_categorized <- read_csv("results/tables/chr21_genes_categorized.csv",
                              show_col_types = FALSE)

eqtl_classification <- read_csv("results/tables/final_gene_classification.csv",
                                show_col_types = FALSE)

cat(sprintf("  Loaded %d chr21 genes\n", nrow(chr21_categorized)))

# Merge eQTL info
chr21_full <- chr21_categorized %>%
  left_join(eqtl_classification %>%
              select(EnsemblID, eQTL_found, eQTL_explains_expression),
            by = "EnsemblID")

cat(sprintf("  Merged eQTL data\n"))

# =============================================================================
# STEP 2: Prepare 5-level alluvial flow data
# =============================================================================

cat("\nStep 2: Preparing alluvial flow data with effect size and eQTL split...\n")

# Effect size threshold for DE genes (based on normalized FC)
# FC < 0.75 = substantial reduction, FC >= 0.75 = modest reduction
EFFECT_SIZE_THRESHOLD <- log2(0.75)  # -0.415

# Create 5-level flow structure with consistent vertical ordering
alluvial_data <- chr21_full %>%
  mutate(
    # Level 1: Starting point
    level1 = "Chr21 Genes",

    # Level 2: FC split (based on raw log2FC - before ploidy correction)
    level2 = ifelse(raw_log2FC >= log2(1.5), ">=1.5 FC", "<1.5 FC"),

    # Level 3: Category split
    # Use actual category names from 02_categorize_genes.R
    level3 = case_when(
      raw_log2FC >= log2(1.5) ~ "Expected\nDosage",
      category == "High Genomic Repeats" ~ "High\nRepeats",
      category == "Low Gene Expression" ~ "Low\nExpression",
      category == "Expression/Ploidy Not Differentially Expressed" ~ "Not DE",
      category == "Expression/Ploidy Differentially Expressed" ~ "DE",
      TRUE ~ "Other"
    ),

    # Level 4: Effect size split for DE genes (based on normalized FC)
    # FC < 0.75 = substantial, FC >= 0.75 = modest
    level4 = case_when(
      level3 != "DE" ~ level3,  # Pass through non-DE categories
      norm_log2FC < EFFECT_SIZE_THRESHOLD ~ "DE: FC<0.75\n(Substantial)",
      TRUE ~ "DE: FC>=0.75\n(Modest)"
    ),

    # Level 5: eQTL split (only for DE genes)
    level5 = case_when(
      level3 != "DE" ~ level4,  # Pass through non-DE categories
      is.na(eQTL_explains_expression) ~ paste0(level4, "\n(No eQTL data)"),
      eQTL_explains_expression == "eQTL supports low expression" ~
        paste0(level4, "\neQTL Supported"),
      eQTL_explains_expression == "No eQTL found" ~
        paste0(level4, "\nNo eQTL"),
      TRUE ~ paste0(level4, "\neQTL Wrong Dir")
    )
  ) %>%
  # Convert to ordered factors to maintain horizontal alignment
  # Order: Expected Dosage at top, then technical categories, then DE categories at bottom
  mutate(
    level2 = factor(level2, levels = c(">=1.5 FC", "<1.5 FC")),
    level3 = factor(level3, levels = c(
      "Expected\nDosage",
      "High\nRepeats",
      "Low\nExpression",
      "Not DE",
      "DE"
    )),
    level4 = factor(level4, levels = c(
      "Expected\nDosage",
      "High\nRepeats",
      "Low\nExpression",
      "Not DE",
      "DE: FC>=0.75\n(Modest)",
      "DE: FC<0.75\n(Substantial)"
    )),
    level5 = factor(level5, levels = c(
      "Expected\nDosage",
      "High\nRepeats",
      "Low\nExpression",
      "Not DE",
      "DE: FC>=0.75\n(Modest)\neQTL Supported",
      "DE: FC>=0.75\n(Modest)\nNo eQTL",
      "DE: FC>=0.75\n(Modest)\neQTL Wrong Dir",
      "DE: FC>=0.75\n(Modest)\n(No eQTL data)",
      "DE: FC<0.75\n(Substantial)\neQTL Supported",
      "DE: FC<0.75\n(Substantial)\nNo eQTL",
      "DE: FC<0.75\n(Substantial)\neQTL Wrong Dir",
      "DE: FC<0.75\n(Substantial)\n(No eQTL data)"
    ))
  ) %>%
  # Count genes in each flow path
  count(level1, level2, level3, level4, level5, name = "n_genes")

cat("  Flow structure (first 10 rows):\n")
print(as.data.frame(head(alluvial_data, 10)))

# Calculate totals
total_genes <- sum(alluvial_data$n_genes)
fc_above <- sum(alluvial_data$n_genes[alluvial_data$level2 == ">=1.5 FC"])
fc_below <- sum(alluvial_data$n_genes[alluvial_data$level2 == "<1.5 FC"])

cat(sprintf("\n  Total chr21 genes: %d\n", total_genes))
cat(sprintf("  FC >= 1.5: %d (%.1f%%)\n",
            fc_above, 100 * fc_above / total_genes))
cat(sprintf("  FC < 1.5: %d (%.1f%%)\n",
            fc_below, 100 * fc_below / total_genes))

# =============================================================================
# STEP 3: Create Panel D using ggalluvial with proper flows
# =============================================================================

cat("\nStep 3: Creating Panel D alluvial plot with ggalluvial...\n")

# Define color palette for final categories (level5)
# Colors are assigned based on the FINAL destination to create progressive effect
color_palette <- c(
  # Expected dosage - cyan/teal
  "Expected\nDosage" = "#5AB4AC",
  # Technical issues - salmon and purple
  "High\nRepeats" = "#E8967A",
  "Low\nExpression" = "#DDA0DD",
  # Not DE - tan
  "Not DE" = "#D2B48C",
  # DE with modest effect (FC >= 0.75) - light reds
  "DE: FC>=0.75\n(Modest)\neQTL Supported" = "#FF9999",
  "DE: FC>=0.75\n(Modest)\nNo eQTL" = "#FF6666",
  "DE: FC>=0.75\n(Modest)\neQTL Wrong Dir" = "#FF3333",
  "DE: FC>=0.75\n(Modest)\n(No eQTL data)" = "#FFCCCC",
  # DE with substantial effect (FC < 0.75) - distinct colors for eQTL status
  "DE: FC<0.75\n(Substantial)\neQTL Supported" = "#4DAF4A",  # Green
  "DE: FC<0.75\n(Substantial)\nNo eQTL" = "#999999",         # Gray
  "DE: FC<0.75\n(Substantial)\neQTL Wrong Dir" = "#FF7F0E",  # Orange
  "DE: FC<0.75\n(Substantial)\n(No eQTL data)" = "#CCCCCC"   # Light gray
)

# Convert to lodes (long) format for proper ggalluvial rendering
alluvial_lodes <- alluvial_data %>%
  mutate(alluvium_id = row_number()) %>%
  to_lodes_form(axes = 1:5, id = "alluvium")

# Create a lookup table for final category (level5) by alluvium
alluvium_to_final <- alluvial_data %>%
  mutate(alluvium = row_number()) %>%
  select(alluvium, final_category = level5)

# Join to get final_category for each alluvium
alluvial_lodes <- alluvial_lodes %>%
  left_join(alluvium_to_final, by = "alluvium")

cat("  Converted to lodes format with final category:\n")
print(head(alluvial_lodes, 15))

# Full color palette including intermediate strata
full_color_palette <- c(
  # Starting point
  "Chr21 Genes" = "#808080",
  # FC split
  ">=1.5 FC" = "#5AB4AC",
  "<1.5 FC" = "#808080",
  # Categories
  "Expected\nDosage" = "#5AB4AC",
  "High\nRepeats" = "#E8967A",
  "Low\nExpression" = "#DDA0DD",
  "Not DE" = "#D2B48C",
  "DE" = "#E41A1C",
  # Effect size
  "DE: FC>=0.75\n(Modest)" = "#FF6B6B",
  "DE: FC<0.75\n(Substantial)" = "#CC0000",
  # Final eQTL categories
  "DE: FC>=0.75\n(Modest)\neQTL Supported" = "#FF9999",
  "DE: FC>=0.75\n(Modest)\nNo eQTL" = "#FF6666",
  "DE: FC>=0.75\n(Modest)\neQTL Wrong Dir" = "#FF3333",
  "DE: FC>=0.75\n(Modest)\n(No eQTL data)" = "#FFCCCC",
  "DE: FC<0.75\n(Substantial)\neQTL Supported" = "#4DAF4A",
  "DE: FC<0.75\n(Substantial)\nNo eQTL" = "#999999",
  "DE: FC<0.75\n(Substantial)\neQTL Wrong Dir" = "#FF7F0E",
  "DE: FC<0.75\n(Substantial)\n(No eQTL data)" = "#CCCCCC"
)

# Create the alluvial plot - flows colored by FINAL destination
p_alluvial <- ggplot(alluvial_lodes,
                     aes(x = x,
                         stratum = stratum,
                         alluvium = alluvium,
                         y = n_genes)) +
  # Draw the flows (alluvia) - colored by final destination
  # Use cubic curve with many segments for smooth S-curves
  # Narrow width to leave room for curves to show
  geom_flow(aes(fill = final_category),
            alpha = 0.7,
            curve_type = "cubic",
            segments = 100,
            width = 1/20) +
  # Draw the strata (narrow bars) - colored by stratum name
  geom_stratum(aes(fill = stratum),
               width = 1/20, color = "black", linewidth = 0.5) +
  # Add text labels
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.0,
            lineheight = 0.8) +
  # Apply color palette

  scale_fill_manual(values = full_color_palette) +
  # X-axis labels
  scale_x_discrete(limits = c("level1", "level2", "level3", "level4", "level5"),
                   labels = c("Chr21\nGenes", "FC\nSplit", "Category",
                              "Effect Size\n(Ploidy-Corrected)", "eQTL\n(FC<0.75)"),
                   expand = c(0.05, 0.05)) +
  labs(
    title = "Chr21 Gene Classification",
    subtitle = "Flows colored by final destination (eQTL category)",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

# Save alluvial plot
ggsave("results/figures/panel_D_with_eqtl.pdf", p_alluvial,
       width = 14, height = 10, bg = "white")
cat("  Saved: results/figures/panel_D_with_eqtl.pdf\n")

# Also save as PNG for quick viewing
ggsave("results/figures/panel_D_with_eqtl.png", p_alluvial,
       width = 14, height = 10, dpi = 150, bg = "white")
cat("  Saved: results/figures/panel_D_with_eqtl.png\n")

# =============================================================================
# STEP 4: Create eQTL breakdown bar plot by effect size
# =============================================================================

cat("\nStep 4: Creating eQTL breakdown bar plot...\n")

# Get breakdown by effect size category
effect_size_breakdown <- alluvial_data %>%
  filter(level3 == "DE") %>%
  group_by(level4) %>%
  summarize(n_genes = sum(n_genes), .groups = "drop") %>%
  mutate(
    percentage = 100 * n_genes / sum(n_genes),
    level4_clean = gsub("\n", " ", level4)
  )

cat("\n  Effect size breakdown:\n")
print(as.data.frame(effect_size_breakdown))

# Create bar plot showing effect size split
p_effect_bar <- ggplot(effect_size_breakdown,
                       aes(x = level4_clean,
                           y = n_genes,
                           fill = level4_clean)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n_genes, percentage)),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "DE: FC<0.75 (Substantial)" = "#D62728",
    "DE: FC>=0.75 (Modest)" = "#FF9896"
  )) +
  labs(
    title = "DE Chr21 Genes by Effect Size",
    subtitle = paste0("Total DE genes: ", sum(effect_size_breakdown$n_genes),
                      " | FC threshold: 0.75"),
    x = "Effect Size Category",
    y = "Number of Genes"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  ylim(0, max(effect_size_breakdown$n_genes) * 1.15)

ggsave("results/figures/eqtl_breakdown_bar.pdf", p_effect_bar,
       width = 8, height = 6)
cat("  Saved: results/figures/eqtl_breakdown_bar.pdf\n")

# =============================================================================
# STEP 5: Create summary table
# =============================================================================

cat("\nStep 5: Creating summary table...\n")

# Overall summary by final category (level5)
summary_table <- alluvial_data %>%
  group_by(level5) %>%
  summarize(count = sum(n_genes), .groups = "drop") %>%
  mutate(
    percentage = 100 * count / sum(count),
    category = gsub("\n", " ", level5)
  ) %>%
  select(category, count, percentage) %>%
  arrange(desc(count))

cat("\nFull category breakdown:\n")
print(as.data.frame(summary_table))

# Save summary
write_csv(summary_table, "results/tables/panel_D_eqtl_summary.csv")
cat("\nSaved: results/tables/panel_D_eqtl_summary.csv\n")

# =============================================================================
# STEP 6: Detailed effect size and eQTL findings summary
# =============================================================================

cat("\n=== Effect Size and eQTL Findings Summary ===\n")

# Count DE genes by effect size
de_genes_total <- sum(alluvial_data$n_genes[alluvial_data$level3 == "DE"])
substantial_effect <- sum(alluvial_data$n_genes[grepl("FC<0.75", alluvial_data$level4)])
modest_effect <- sum(alluvial_data$n_genes[grepl("FC>=0.75", alluvial_data$level4)])

cat(sprintf("\nDE genes (significant after ploidy correction): %d\n", de_genes_total))
cat(sprintf("  Substantial effect (FC < 0.75): %d (%.1f%%)\n",
            substantial_effect, 100 * substantial_effect / de_genes_total))
cat(sprintf("  Modest effect (FC >= 0.75):     %d (%.1f%%)\n",
            modest_effect, 100 * modest_effect / de_genes_total))

# eQTL breakdown within each effect size category
cat("\n  eQTL support by effect size:\n")

eqtl_by_effect <- alluvial_data %>%
  filter(level3 == "DE") %>%
  mutate(
    effect_category = ifelse(grepl("FC<0.75", level4), "Substantial", "Modest"),
    eqtl_status = case_when(
      grepl("eQTL Supported", level5) ~ "eQTL Supported",
      grepl("No eQTL", level5) ~ "No eQTL",
      grepl("Wrong Dir", level5) ~ "eQTL Wrong Direction",
      TRUE ~ "No eQTL data"
    )
  ) %>%
  group_by(effect_category, eqtl_status) %>%
  summarize(n = sum(n_genes), .groups = "drop")

print(as.data.frame(eqtl_by_effect))

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat(sprintf("\nTotal chr21 protein-coding genes: %d\n", total_genes))
cat(sprintf("Expected Dosage (raw FC >= 1.5):  %d (%.1f%%)\n",
            fc_above, 100 * fc_above / total_genes))
cat(sprintf("Below Expected (raw FC < 1.5):    %d (%.1f%%)\n",
            fc_below, 100 * fc_below / total_genes))

# Breakdown of <1.5 FC genes
not_de <- sum(alluvial_data$n_genes[alluvial_data$level3 == "Not DE"])
low_expr <- sum(alluvial_data$n_genes[alluvial_data$level3 == "Low\nExpression"])
high_rep <- sum(alluvial_data$n_genes[alluvial_data$level3 == "High\nRepeats"])

cat(sprintf("\nBreakdown of FC < 1.5 genes:\n"))
cat(sprintf("  Not DE (padj >= 0.05):      %d\n", not_de))
cat(sprintf("  Low Expression:             %d\n", low_expr))
cat(sprintf("  High Repeats:               %d\n", high_rep))
cat(sprintf("  DE (padj < 0.05):           %d\n", de_genes_total))
cat(sprintf("    - Substantial (FC<0.75):  %d\n", substantial_effect))
cat(sprintf("    - Modest (FC>=0.75):      %d\n", modest_effect))

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/figures/alluvial_eqtl_session_info.txt")

cat("\n=== Enhanced Alluvial Plot Complete ===\n")
cat("Figures saved to results/figures/\n")
cat("  - panel_D_with_eqtl.pdf (5-level alluvial with effect size)\n")
cat("  - panel_D_with_eqtl.png (PNG version)\n")
cat("  - eqtl_breakdown_bar.pdf (effect size breakdown)\n\n")
