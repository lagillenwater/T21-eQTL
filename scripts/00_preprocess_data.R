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
# STEP 3b: Analysis-cohort definition and Table 1
# =============================================================================
#
# The cohort is asymmetric: T21 needs RNA-seq AND chr21 WGS (genotypes feed
# the within-T21 dosage regressions only); Control needs RNA-seq alone. See
# scripts/lib/cohort.R for the full rationale.

cat("\nStep 3b: Defining analysis cohort and building Table 1...\n")

source("scripts/lib/cohort.R")
source("scripts/lib/table1.R")

if (!dir.exists("results/tables")) {
  dir.create("results/tables", recursive = TRUE)
}

# WGS roster from the PASS file headers only - cheap, and needed here because
# script 01 must run DE on the analysis cohort, before script 02 processes any
# genotypes.
wgs <- wgs_subjects(c("data/chr21_ds_PASS.csv", "data/chr21_ctrl_PASS.csv"))
cat(sprintf("  Subjects with chr21 WGS (from headers): %d\n", length(wgs)))

metadata_matched <- metadata_matched %>%
  mutate(subject_id = subject_id_from_labid(LabID),
         has_wgs = subject_id %in% wgs)

print(as.data.frame(metadata_matched %>%
  count(Karyotype, name = "rnaseq") %>%
  left_join(metadata_matched %>% filter(has_wgs) %>% count(Karyotype, name = "with_wgs"),
            by = "Karyotype")))

cohort <- analysis_cohort(as.data.table(metadata_matched))
cat(sprintf("  Analysis cohort: T21 %d, Control %d, total %d\n",
            sum(cohort$Karyotype == "T21"), sum(cohort$Karyotype == "Control"),
            nrow(cohort)))
stopifnot(sum(cohort$Karyotype == "T21") == 302,
          sum(cohort$Karyotype == "Control") == 95)
fwrite(cohort, "data/processed/analysis_cohort.csv")
if (!file.exists("data/processed/analysis_cohort.csv")) {
  stop("failed to write data/processed/analysis_cohort.csv")
}
cat("  Saved: data/processed/analysis_cohort.csv\n")

# --- Table 1 -----------------------------------------------------------

t21a <- cohort[Karyotype == "T21"]
ctla <- cohort[Karyotype == "Control"]

row_cont <- function(label, var, digits = 1) {
  cmp <- compare_groups(cohort[[var]], cohort$Karyotype)
  data.table(characteristic = label, level = "median [IQR]",
             t21 = summarize_continuous(t21a[[var]], digits),
             control = summarize_continuous(ctla[[var]], digits),
             overall = summarize_continuous(cohort[[var]], digits),
             p_value = cmp$p, test = cmp$test,
             n_missing = sum(is.na(cohort[[var]])))
}
row_cat <- function(label, var) {
  cmp  <- compare_groups(cohort[[var]], cohort$Karyotype)
  levs <- sort(unique(as.character(cohort[[var]][!is.na(cohort[[var]])])))
  rbindlist(lapply(seq_along(levs), function(i)
    data.table(characteristic = if (i == 1) label else "", level = levs[i],
               t21 = summarize_categorical(t21a[[var]], levs[i]),
               control = summarize_categorical(ctla[[var]], levs[i]),
               overall = summarize_categorical(cohort[[var]], levs[i]),
               p_value = if (i == 1) cmp$p else NA_real_,
               test = if (i == 1) cmp$test else "",
               n_missing = if (i == 1) sum(is.na(cohort[[var]])) else NA_integer_)))
}

tbl1 <- rbindlist(list(
  data.table(characteristic = "Subjects in analysis cohort", level = "n",
             t21 = as.character(nrow(t21a)), control = as.character(nrow(ctla)),
             overall = as.character(nrow(cohort)),
             p_value = NA_real_, test = "", n_missing = NA_integer_),
  data.table(characteristic = "Excluded: RNA-seq without WGS", level = "n",
             t21 = as.character(sum(metadata_matched$Karyotype == "T21" & !metadata_matched$has_wgs)),
             control = "0 (not required)",
             overall = as.character(sum(metadata_matched$Karyotype == "T21" & !metadata_matched$has_wgs)),
             p_value = NA_real_, test = "", n_missing = NA_integer_),
  row_cont("Age at visit (years)", "Age_at_visit"),
  row_cont("BMI", "BMI"),
  row_cat("Sex", "Sex"),
  row_cat("Sample source", "Sample_source"),
  row_cat("Study visit", "Event_name")))

# Optional comorbidity block. Long format with ONE comment line before the
# header: RecordID, Condition, HasCondition, Age.group, min_Age, max_Age.
como_path <- "data/P4C_Comorbidity_020921.tsv"
if (file.exists(como_path)) {
  como <- fread(como_path, skip = 1)
  if (all(c("RecordID", "Condition", "HasCondition") %in% names(como))) {
    como <- como[RecordID %in% cohort$RecordID]
    top  <- como[, .(n_with = sum(HasCondition == 1, na.rm = TRUE)), by = Condition][
      order(-n_with)][seq_len(min(5, .N))]
    wide <- dcast(como, RecordID ~ Condition, value.var = "HasCondition")
    ann  <- merge(cohort[, .(RecordID, Karyotype)], wide, by = "RecordID")
    tbl1 <- rbindlist(list(tbl1, rbindlist(lapply(top$Condition, function(cn) {
      cmp <- compare_groups(as.character(ann[[cn]]), ann$Karyotype)
      data.table(characteristic = cn, level = "n (%) with condition",
                 t21 = summarize_categorical(ann[Karyotype == "T21"][[cn]], 1),
                 control = summarize_categorical(ann[Karyotype == "Control"][[cn]], 1),
                 overall = summarize_categorical(ann[[cn]], 1),
                 p_value = cmp$p, test = cmp$test, n_missing = sum(is.na(ann[[cn]])))
    }))))
  } else {
    cat("  comorbidity file present but unexpected columns; skipping that block\n")
  }
}

tbl1[, p_value := ifelse(is.na(p_value), "", format.pval(p_value, digits = 2, eps = 1e-4))]
fwrite(tbl1, "results/tables/table1_analysis_cohort.csv")
writeLines(c(
  "# Table 1. Characteristics of the analysis cohort", "",
  sprintf("T21 subjects require both whole-blood RNA-seq and chr21 WGS (n = %d of %d).",
          nrow(t21a), sum(metadata_matched$Karyotype == "T21")),
  sprintf("Controls require RNA-seq only (n = %d), because genotypes are used solely",
          nrow(ctla)),
  "for the within-T21 dosage regressions, which controls do not enter.",
  "Race and ethnicity are not recorded in the available metadata.", "",
  paste("|", paste(names(tbl1), collapse = " | "), "|"),
  paste("|", paste(rep("---", ncol(tbl1)), collapse = " | "), "|"),
  apply(tbl1, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))),
  "results/tables/table1_analysis_cohort.md")
for (f in c("results/tables/table1_analysis_cohort.csv",
            "results/tables/table1_analysis_cohort.md")) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat("  Wrote ", f, "\n", sep = "")
}
print(tbl1)

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

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  ADDED the analysis-cohort definition and Table 1.
#             The cohort is asymmetric by design: T21 require RNA-seq AND chr21
#             WGS (302 of 304), Controls require RNA-seq only (95 of 95).
#             Reason: genotypes are used only for the within-T21 dosage
#             regressions, which controls never enter, so requiring WGS of a
#             control would discard 89 of 95 for no analytic gain. Dropping the
#             2 ungenotyped T21 makes the DE cohort and the eQTL cohort the same
#             people; previously DE ran on 304 and eQTL on 302.
#             The WGS roster is read from the PASS file HEADERS only, so this
#             stays cheap enough to run before any genotype processing.
#             sample_metadata.csv keeps all 399 rows and gains subject_id and
#             has_wgs; data/processed/analysis_cohort.csv carries the 397-row
#             roster that scripts 01-04 subset to.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
