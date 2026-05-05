# 04_alluvial_plot.R
#
# Purpose: Create alluvial/Sankey diagram reproducing Panel D from Hunter et al
#          Shows flow of chr21 genes through categorization
#
# Inputs:
#   - results/tables/chr21_genes_categorized.csv
#   - results/tables/category_summary.csv
#
# Outputs:
#   - results/figures/panel_D_alluvial.pdf
#   - results/figures/panel_D_barplot.pdf (alternative visualization)
#
# Author: Claude Code
# Date: 2025-11-11

# Load required libraries
library(tidyverse)
library(ggalluvial)

# Set seed for reproducibility
set.seed(42)

cat("=== T21-eQTL Analysis: Alluvial Plot (Panel D) ===\n\n")

# =============================================================================
# STEP 1: Load categorized data
# =============================================================================

cat("Step 1: Loading categorized chr21 genes...\n")

if (!file.exists("results/tables/chr21_genes_categorized.csv")) {
  stop("Categorized data not found. Run 02_categorize_genes.R first.")
}

chr21_categorized <- read_csv("results/tables/chr21_genes_categorized.csv",
                              show_col_types = FALSE)

cat(sprintf("  Loaded %d chr21 genes\n", nrow(chr21_categorized)))

# =============================================================================
# STEP 2: Prepare alluvial flow data
# =============================================================================

cat("\nStep 2: Preparing alluvial flow data...\n")

# Create 3-level flow structure matching Panel D
alluvial_data <- chr21_categorized %>%
  mutate(
    # Level 1: Starting point
    level1 = "Chr21 Genes",

    # Level 2: FC split
    level2 = ifelse(raw_log2FC >= log2(1.5), ">=1.5 FC", "<1.5 FC"),

    # Level 3: Final categories
    level3 = case_when(
      raw_log2FC >= log2(1.5) ~ ">=1.5 FC",
      category == "High Genomic Repeats" ~ "High Genomic\nRepeats",
      category == "Low Gene Expression" ~ "Low Gene\nExpression",
      category == "Not Differentially Expressed" ~ "Not\nDifferentially Expressed",
      category == "Differentially Expressed (<1.5 FC)" ~ "Differentially\nExpressed",
      TRUE ~ "Other"
    )
  ) %>%
  # Count genes in each flow path
  count(level1, level2, level3, name = "n_genes")

cat("  Flow structure:\n")
print(as.data.frame(alluvial_data))

# Calculate totals for labels
total_genes <- sum(alluvial_data$n_genes)
fc_above <- sum(alluvial_data$n_genes[alluvial_data$level2 == ">=1.5 FC"])
fc_below <- sum(alluvial_data$n_genes[alluvial_data$level2 == "<1.5 FC"])

cat(sprintf("\n  Total chr21 genes: %d\n", total_genes))
cat(sprintf("  FC >= 1.5: %d (%.1f%%)\n",
            fc_above, 100 * fc_above / total_genes))
cat(sprintf("  FC < 1.5: %d (%.1f%%)\n",
            fc_below, 100 * fc_below / total_genes))

# Print detailed breakdown
cat("\n  Final category breakdown:\n")
category_counts <- alluvial_data %>%
  group_by(level3) %>%
  summarize(n = sum(n_genes), .groups = "drop") %>%
  arrange(desc(n))

for (i in 1:nrow(category_counts)) {
  cat(sprintf("    %s: %d (%.1f%%)\n",
              gsub("\n", " ", category_counts$level3[i]),
              category_counts$n[i],
              100 * category_counts$n[i] / total_genes))
}

# =============================================================================
# STEP 3: Create alluvial plot
# =============================================================================

cat("\nStep 3: Creating alluvial diagram...\n")

# Color scheme matching paper Panel D
category_colors <- c(
  ">=1.5 FC" = "#5AB4AC",  # Cyan/teal - Expected dosage
  "High Genomic\nRepeats" = "#D8B365",  # Tan/beige
  "Low Gene\nExpression" = "#C77CFF",   # Purple/lavender
  "Not\nDifferentially Expressed" = "#AF8DC3",  # Mauve/brown
  "Differentially\nExpressed" = "#E08214"  # Orange
)

# Create alluvial plot matching Panel D style with 3-level structure
p_alluvial <- ggplot(alluvial_data,
                     aes(axis1 = level1,
                         axis2 = level2,
                         axis3 = level3,
                         y = n_genes)) +
  # Add smooth alluvial flows with transparency
  geom_alluvium(aes(fill = level3),
                width = 1/12,
                alpha = 0.75,
                curve_type = "sigmoid") +
  # Add stratum rectangles
  geom_stratum(width = 1/12,
               aes(fill = level3),
               color = "white",
               linewidth = 0.5) +
  # Add text labels on strata
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3,
            fontface = "bold",
            lineheight = 0.85) +
  # Colors matching paper
  scale_fill_manual(values = category_colors,
                    name = "Category",
                    labels = gsub("\n", " ", names(category_colors))) +
  scale_x_discrete(limits = c("Total", "FC Split", "Categories"),
                   expand = c(0.05, 0.05)) +
  labs(
    title = "RNA-seq Accounting For Trisomy",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.position = "none",  # Remove legend, labels are in plot
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save alluvial plot
ggsave("results/figures/panel_D_alluvial.pdf", p_alluvial,
       width = 10, height = 7, bg = "white")
cat("  Saved: results/figures/panel_D_alluvial.pdf\n")

# =============================================================================
# STEP 4: Create alternative bar plot visualization
# =============================================================================

cat("\nStep 4: Creating alternative bar plot...\n")

# Prepare data for bar plot
bar_data <- alluvial_data %>%
  group_by(level3) %>%
  summarize(n_genes = sum(n_genes), .groups = "drop") %>%
  mutate(
    # Create ordered factor for better visualization
    level3_clean = gsub("\n", " ", level3),
    level3_clean = factor(level3_clean,
                          levels = c(">=1.5 FC",
                                     "High Genomic Repeats",
                                     "Low Gene Expression",
                                     "Not Differentially Expressed",
                                     "Differentially Expressed"))
  )

# Create bar plot
p_bar <- ggplot(bar_data, aes(x = level3_clean, y = n_genes,
                               fill = level3)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_genes),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = category_colors,
                    name = "Category") +
  labs(
    title = "Chr21 Gene Categories",
    subtitle = paste0("Total genes analyzed: ", total_genes),
    x = "Category",
    y = "Number of Genes"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )

# Save bar plot
ggsave("results/figures/panel_D_barplot.pdf", p_bar,
       width = 10, height = 6)
cat("  Saved: results/figures/panel_D_barplot.pdf\n")

# =============================================================================
# STEP 5: Create detailed category breakdown
# =============================================================================

cat("\nStep 5: Creating detailed category breakdown plot...\n")

# Calculate percentages for final categories
category_breakdown <- alluvial_data %>%
  group_by(level3) %>%
  summarize(n_genes = sum(n_genes), .groups = "drop") %>%
  mutate(
    percentage = 100 * n_genes / sum(n_genes),
    level3_clean = gsub("\n", " ", level3)
  ) %>%
  arrange(desc(n_genes))

# Create pie/donut chart
p_breakdown <- ggplot(category_breakdown,
                      aes(x = "", y = n_genes, fill = level3)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n_genes, percentage)),
            position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = category_colors,
                    name = "Category") +
  labs(
    title = "Chr21 Gene Category Distribution",
    subtitle = paste0("Total: ", total_genes, " genes")
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm")
  )

ggsave("results/figures/panel_D_pie_chart.pdf", p_breakdown,
       width = 10, height = 6)
cat("  Saved: results/figures/panel_D_pie_chart.pdf\n")

# =============================================================================
# STEP 6: Summary table
# =============================================================================

cat("\n=== Panel D Summary ===\n")

summary_table <- category_breakdown %>%
  select(Category = level3_clean, Count = n_genes, Percentage = percentage) %>%
  arrange(desc(Count))

cat("\nCategory breakdown:\n")
print(as.data.frame(summary_table))

# Compare with paper
cat("\n=== Comparison with Hunter et al. (2023) ===\n")
cat("Expected proportions from paper:\n")
cat("  Expected Dosage: ~38.5% (101/262)\n")
cat("  High Genomic Repeats: ~3.4% (9/262)\n")
cat("  Low Expression: ~22.1% (58/262)\n")
cat("  Not DE: ~30.2% (79/262)\n")
cat("  DE (<1.5 FC): ~1.9% (5/262)\n")

cat("\n=== Key Finding ===\n")
expected_pct <- summary_table %>%
  filter(grepl(">=1.5 FC", Category)) %>%
  pull(Percentage) %>%
  head(1)  # Take first match

de_count <- summary_table %>%
  filter(grepl("Differentially Expressed", Category)) %>%
  pull(Count) %>%
  head(1)  # Take first match

if (length(expected_pct) > 0 && !is.na(expected_pct) && expected_pct > 30) {
  cat(sprintf("%.1f%% of chr21 genes show expected 1.5-fold change.\n",
              expected_pct))
  if (expected_pct < 50) {
    cat("Most chr21 genes show lower than expected expression.\n")
  } else {
    cat("This provides evidence against widespread dosage compensation.\n")
  }
}

if (length(de_count) > 0 && !is.na(de_count) && de_count <= 10) {
  cat(sprintf("\nOnly %d genes show significantly lower expression.\n",
              de_count))
  cat("These will be investigated for eQTL explanations.\n")
} else if (length(de_count) > 0 && !is.na(de_count)) {
  cat(sprintf("\n%d genes show significantly lower expression.\n", de_count))
}

# Save summary table
write_csv(summary_table, "results/tables/panel_D_summary.csv")
cat("\nSaved: results/tables/panel_D_summary.csv\n")

# Save session info
writeLines(capture.output(sessionInfo()),
           "results/figures/alluvial_plot_session_info.txt")

cat("\n=== Alluvial Plot Complete ===\n")
cat("Next step: Run 05_eqtl_analysis.R\n\n")
