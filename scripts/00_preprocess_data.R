# 00_preprocess_data.R
#
# Purpose: Convert long-format RNA-seq count data to gene x sample matrix
#          Match samples with metadata and create gene annotation table
#
# Inputs:
#   - data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt (long format)
#   - data/P4C_metadata_021921_Costello.txt
#
# Outputs:
#   - data/processed/count_matrix.csv
#   - data/processed/sample_metadata.csv
#   - data/processed/gene_annotations.csv
#
# Date: 2025-11-11

# Load required libraries
library(tidyverse)
library(data.table)

# Set seed for reproducibility
set.seed(42)

cat("=== T21-eQTL Analysis: Data Preprocessing ===\n\n")

# Create output directory if needed
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
  cat("Created data/processed/ directory\n")
}

# =============================================================================
# STEP 1: Load and inspect raw count data
# =============================================================================

cat("Step 1: Loading raw count data...\n")

# Check that input file exists
if (!file.exists("data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt")) {
  stop("Count data file not found: data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt")
}

# Load count data using data.table for efficiency (large file)
count_data <- fread("data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt",
                    select = c("LabID", "EnsemblID", "Gene_name",
                               "Chr", "Gene_type", "Value"))

cat(sprintf("  Loaded %d rows\n", nrow(count_data)))
cat(sprintf("  Unique samples: %d\n", length(unique(count_data$LabID))))
cat(sprintf("  Unique genes: %d\n", length(unique(count_data$EnsemblID))))

# =============================================================================
# STEP 2: Load and process metadata
# =============================================================================

cat("\nStep 2: Loading sample metadata...\n")

if (!file.exists("data/P4C_metadata_021921_Costello.txt")) {
  stop("Metadata file not found: data/P4C_metadata_021921_Costello.txt")
}

metadata <- read_tsv("data/P4C_metadata_021921_Costello.txt",
                     show_col_types = FALSE)

cat(sprintf("  Loaded %d samples\n", nrow(metadata)))

# Check karyotype distribution
karyotype_counts <- table(metadata$Karyotype)
cat("  Karyotype distribution:\n")
print(karyotype_counts)

# Filter to only T21 and Control
metadata <- metadata %>%
  filter(Karyotype %in% c("T21", "Control"))

cat(sprintf("  Kept %d samples (T21: %d, Control: %d)\n",
            nrow(metadata),
            sum(metadata$Karyotype == "T21"),
            sum(metadata$Karyotype == "Control")))

# =============================================================================
# STEP 3: Match sample IDs between counts and metadata
# =============================================================================

cat("\nStep 3: Matching sample IDs...\n")

# Check which samples have both count and metadata
# The LabIDs in both files should match directly (no stripping needed)
samples_in_counts <- unique(count_data$LabID)
samples_in_metadata <- metadata$LabID

samples_both <- intersect(samples_in_counts, samples_in_metadata)
samples_counts_only <- setdiff(samples_in_counts, samples_in_metadata)
samples_metadata_only <- setdiff(samples_in_metadata, samples_in_counts)

cat(sprintf("  Samples in both count and metadata: %d\n",
            length(samples_both)))
cat(sprintf("  Samples in counts only: %d\n",
            length(samples_counts_only)))
cat(sprintf("  Samples in metadata only: %d\n",
            length(samples_metadata_only)))

# Show examples of unmatched samples
if (length(samples_counts_only) > 0) {
  cat("  Example count-only samples:\n")
  print(head(samples_counts_only, 5))
}
if (length(samples_metadata_only) > 0) {
  cat("  Example metadata-only samples:\n")
  print(head(samples_metadata_only, 5))
}

# Filter counts to only samples with metadata
count_data_matched <- count_data %>%
  filter(LabID %in% samples_both)

# Filter metadata to only samples with counts
metadata_matched <- metadata %>%
  filter(LabID %in% samples_both)

cat(sprintf("  Final matched samples: %d\n", length(samples_both)))

# =============================================================================
# STEP 4: Create gene annotation table
# =============================================================================

cat("\nStep 4: Creating gene annotation table...\n")

gene_annotations <- count_data_matched %>%
  select(EnsemblID, Gene_name, Chr, Gene_type) %>%
  distinct()

cat(sprintf("  Total unique genes: %d\n", nrow(gene_annotations)))

# Count genes by chromosome
chr_counts <- gene_annotations %>%
  count(Chr) %>%
  arrange(desc(n))

cat("  Genes per chromosome (top 10):\n")
print(head(chr_counts, 10))

# Specifically check chr21
chr21_count <- sum(gene_annotations$Chr == "chr21")
cat(sprintf("  Chr21 genes: %d\n", chr21_count))

# =============================================================================
# STEP 5: Convert to wide format (gene x sample matrix)
# =============================================================================

cat("\nStep 5: Converting to wide format matrix...\n")
cat("  This may take several minutes for large datasets...\n")

# Pivot to wide format: rows = genes, columns = samples
count_matrix <- count_data_matched %>%
  select(EnsemblID, LabID, Value) %>%
  pivot_wider(names_from = LabID,
              values_from = Value,
              values_fill = 0) %>%
  column_to_rownames("EnsemblID")

cat(sprintf("  Count matrix dimensions: %d genes x %d samples\n",
            nrow(count_matrix), ncol(count_matrix)))

# =============================================================================
# STEP 6: Quality control checks
# =============================================================================

cat("\nStep 6: Quality control checks...\n")

# Check for NA values
na_count <- sum(is.na(count_matrix))
if (na_count > 0) {
  warning(sprintf("Count matrix contains %d NA values", na_count))
} else {
  cat("  No NA values in count matrix\n")
}

# Check that all samples in matrix are in metadata
matrix_samples <- colnames(count_matrix)
metadata_samples <- metadata_matched$LabID
unmatched <- setdiff(matrix_samples, metadata_samples)
if (length(unmatched) > 0) {
  warning(sprintf("%d samples in matrix not in metadata",
                  length(unmatched)))
} else {
  cat("  All matrix samples found in metadata\n")
}

# Check gene filtering - remove very lowly expressed genes
# Calculate mean expression per gene
gene_means <- rowMeans(count_matrix)
low_expr_genes <- sum(gene_means < 1)
cat(sprintf("  Genes with mean count < 1: %d (%.1f%%)\n",
            low_expr_genes, 100 * low_expr_genes / nrow(count_matrix)))

# Filter genes with mean count >= 1
genes_to_keep <- gene_means >= 1
count_matrix_filtered <- count_matrix[genes_to_keep, ]
gene_annotations_filtered <- gene_annotations %>%
  filter(EnsemblID %in% rownames(count_matrix_filtered))

cat(sprintf("  Filtered count matrix: %d genes x %d samples\n",
            nrow(count_matrix_filtered), ncol(count_matrix_filtered)))

# Check chr21 genes after filtering
chr21_after_filter <- sum(gene_annotations_filtered$Chr == "chr21")
cat(sprintf("  Chr21 genes after filtering: %d\n", chr21_after_filter))

# =============================================================================
# STEP 7: Prepare and save outputs
# =============================================================================

cat("\nStep 7: Saving processed data...\n")

# Add gene info to count matrix for easier reading
count_matrix_with_info <- count_matrix_filtered %>%
  rownames_to_column("EnsemblID") %>%
  left_join(gene_annotations_filtered %>%
              select(EnsemblID, Gene_name, Chr),
            by = "EnsemblID") %>%
  select(EnsemblID, Gene_name, Chr, everything())

# Save count matrix
write_csv(count_matrix_with_info,
          "data/processed/count_matrix.csv")
cat("  Saved: data/processed/count_matrix.csv\n")

# Save metadata
metadata_matched %>%
  write_csv("data/processed/sample_metadata.csv")
cat("  Saved: data/processed/sample_metadata.csv\n")

# Save gene annotations
gene_annotations_filtered %>%
  write_csv("data/processed/gene_annotations.csv")
cat("  Saved: data/processed/gene_annotations.csv\n")

# =============================================================================
# STEP 8: Summary statistics
# =============================================================================

cat("\n=== Preprocessing Summary ===\n")
cat(sprintf("Samples: %d total (%d T21, %d Control)\n",
            nrow(metadata_matched),
            sum(metadata_matched$Karyotype == "T21"),
            sum(metadata_matched$Karyotype == "Control")))
cat(sprintf("Genes: %d total, %d on chr21\n",
            nrow(gene_annotations_filtered),
            chr21_after_filter))
cat(sprintf("Count matrix: %d x %d\n",
            nrow(count_matrix_filtered),
            ncol(count_matrix_filtered)))

# Sample breakdown by sex and karyotype
sample_breakdown <- metadata_matched %>%
  count(Karyotype, Sex) %>%
  arrange(Karyotype, Sex)

cat("\nSample breakdown by karyotype and sex:\n")
print(as.data.frame(sample_breakdown))

# Gene type breakdown for chr21
if (chr21_after_filter > 0) {
  chr21_gene_types <- gene_annotations_filtered %>%
    filter(Chr == "chr21") %>%
    count(Gene_type) %>%
    arrange(desc(n))

  cat("\nChr21 gene types:\n")
  print(as.data.frame(chr21_gene_types))
}

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()),
           "data/processed/preprocessing_session_info.txt")
cat("\nSaved session info to data/processed/preprocessing_session_info.txt\n")

cat("\n=== Preprocessing Complete ===\n")
cat("Next step: Run 01_deseq2_analysis.R\n\n")
