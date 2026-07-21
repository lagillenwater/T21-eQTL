# process_blacklist.R
#
# Purpose: Identify chr21 genes that overlap with ENCODE blacklist regions
#          These genes have high genomic repeats causing mapping artifacts
#
# Inputs:
#   - data/hg38-blacklist.v2.bed.gz (ENCODE blacklist)
#   - data/gencode.v44.basic.annotation.gtf.gz (gene coordinates)
#
# Outputs:
#   - data/processed/blacklisted_genes.csv
#
# Date: 2025-12-04

library(tidyverse)
library(GenomicRanges)

cat("=== Processing ENCODE Blacklist for Chr21 Genes ===\n\n")

# =============================================================================
# STEP 1: Load blacklist regions
# =============================================================================

cat("Step 1: Loading blacklist regions...\n")

blacklist_file <- "data/hg38-blacklist.v2.bed.gz"
if (!file.exists(blacklist_file)) {
  stop("Blacklist file not found: ", blacklist_file)
}

blacklist <- read_tsv(blacklist_file,
                      col_names = c("chr", "start", "end", "type"),
                      col_types = "ciic")

cat(sprintf("  Loaded %d blacklist regions\n", nrow(blacklist)))

# Filter to chr21
blacklist_chr21 <- blacklist %>%
  filter(chr == "chr21")

cat(sprintf("  Chr21 blacklist regions: %d\n", nrow(blacklist_chr21)))
cat("\n  Chr21 blacklist regions:\n")
print(as.data.frame(blacklist_chr21))

# Create GRanges object for blacklist
blacklist_gr <- GRanges(
  seqnames = blacklist_chr21$chr,
  ranges = IRanges(start = blacklist_chr21$start + 1,  # BED is 0-based
                   end = blacklist_chr21$end),
  type = blacklist_chr21$type
)

# =============================================================================
# STEP 2: Load gene annotations from GENCODE GTF
# =============================================================================

cat("\nStep 2: Loading gene annotations from GENCODE...\n")

gtf_file <- "data/gencode.v44.basic.annotation.gtf.gz"
if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

# Parse GTF for chr21 genes only (more efficient than loading entire file)
gtf_lines <- read_lines(gtf_file)
gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]  # Remove comments

# Filter to chr21 gene entries
chr21_gene_lines <- gtf_lines[grepl("^chr21\t.*\tgene\t", gtf_lines)]
cat(sprintf("  Found %d chr21 gene entries\n", length(chr21_gene_lines)))

# Parse gene information
parse_gtf_gene <- function(line) {
  fields <- strsplit(line, "\t")[[1]]
  chr <- fields[1]
  start <- as.integer(fields[4])
  end <- as.integer(fields[5])
  strand <- fields[7]
  attributes <- fields[9]

  # Extract gene_id and gene_name from attributes
  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attributes)
  gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attributes)
  gene_type <- sub('.*gene_type "([^"]+)".*', "\\1", attributes)

  # Remove version from gene_id for matching
  gene_id_base <- sub("\\..*", "", gene_id)

  data.frame(
    chr = chr,
    start = start,
    end = end,
    strand = strand,
    gene_id = gene_id,
    gene_id_base = gene_id_base,
    gene_name = gene_name,
    gene_type = gene_type,
    stringsAsFactors = FALSE
  )
}

chr21_genes <- do.call(rbind, lapply(chr21_gene_lines, parse_gtf_gene))
cat(sprintf("  Parsed %d chr21 genes\n", nrow(chr21_genes)))

# Create GRanges for genes
genes_gr <- GRanges(
  seqnames = chr21_genes$chr,
  ranges = IRanges(start = chr21_genes$start, end = chr21_genes$end),
  strand = chr21_genes$strand,
  gene_id = chr21_genes$gene_id,
  gene_id_base = chr21_genes$gene_id_base,
  gene_name = chr21_genes$gene_name,
  gene_type = chr21_genes$gene_type
)

# =============================================================================
# STEP 3: Find overlapping genes
# =============================================================================

cat("\nStep 3: Finding genes overlapping blacklist regions...\n")

# Find overlaps
overlaps <- findOverlaps(genes_gr, blacklist_gr)

# Get blacklisted genes
blacklisted_idx <- unique(queryHits(overlaps))
blacklisted_genes <- chr21_genes[blacklisted_idx, ]

cat(sprintf("  Genes overlapping blacklist: %d\n", nrow(blacklisted_genes)))

# Calculate overlap fraction for each gene
blacklisted_genes$overlap_bp <- 0
blacklisted_genes$overlap_fraction <- 0
blacklisted_genes$blacklist_type <- ""

for (i in seq_len(nrow(blacklisted_genes))) {
  gene_gr <- genes_gr[blacklisted_idx[i]]
  gene_length <- width(gene_gr)

  # Find which blacklist regions overlap this gene
  gene_overlaps <- findOverlaps(gene_gr, blacklist_gr)

  if (length(gene_overlaps) > 0) {
    # Get overlapping blacklist regions
    bl_hits <- subjectHits(gene_overlaps)
    overlapping_bl <- blacklist_gr[bl_hits]

    # Calculate intersection manually
    total_overlap <- 0
    for (j in seq_along(overlapping_bl)) {
      bl_region <- overlapping_bl[j]
      overlap_start <- max(start(gene_gr), start(bl_region))
      overlap_end <- min(end(gene_gr), end(bl_region))
      if (overlap_end >= overlap_start) {
        total_overlap <- total_overlap + (overlap_end - overlap_start + 1)
      }
    }

    blacklisted_genes$overlap_bp[i] <- total_overlap
    blacklisted_genes$overlap_fraction[i] <- total_overlap / gene_length

    # Get blacklist types
    bl_types <- unique(blacklist_gr$type[bl_hits])
    blacklisted_genes$blacklist_type[i] <- paste(bl_types, collapse = ";")
  }
}

# =============================================================================
# STEP 4: Filter to protein-coding genes with significant overlap
# =============================================================================

cat("\nStep 4: Filtering results...\n")

# Keep protein-coding genes with at least 10% overlap
MIN_OVERLAP_FRACTION <- 0.10

blacklisted_protein_coding <- blacklisted_genes %>%
  filter(gene_type == "protein_coding") %>%
  filter(overlap_fraction >= MIN_OVERLAP_FRACTION) %>%
  arrange(desc(overlap_fraction))

cat(sprintf("  Protein-coding genes with >= %.0f%% overlap: %d\n",
            MIN_OVERLAP_FRACTION * 100, nrow(blacklisted_protein_coding)))

if (nrow(blacklisted_protein_coding) > 0) {
  cat("\n  Blacklisted protein-coding genes:\n")
  print(blacklisted_protein_coding %>%
          select(gene_name, gene_id_base, start, end,
                 overlap_bp, overlap_fraction, blacklist_type) %>%
          mutate(overlap_fraction = round(overlap_fraction, 3)) %>%
          as.data.frame())
}

# Also show all blacklisted genes (any type)
cat(sprintf("\n  All gene types with >= %.0f%% overlap: %d\n",
            MIN_OVERLAP_FRACTION * 100,
            sum(blacklisted_genes$overlap_fraction >= MIN_OVERLAP_FRACTION)))

# =============================================================================
# STEP 5: Load existing gene annotations and match
# =============================================================================

cat("\nStep 5: Matching to analysis gene list...\n")

gene_annotations <- read_csv("data/processed/gene_annotations.csv",
                             show_col_types = FALSE)

# Extract base Ensembl ID (without version)
gene_annotations <- gene_annotations %>%
  mutate(gene_id_base = sub("\\..*", "", EnsemblID))

chr21_analysis_genes <- gene_annotations %>%
  filter(Chr == "chr21", Gene_type == "protein_coding")

cat(sprintf("  Chr21 protein-coding genes in analysis: %d\n",
            nrow(chr21_analysis_genes)))

# Match blacklisted genes to analysis genes
blacklisted_in_analysis <- chr21_analysis_genes %>%
  filter(gene_id_base %in% blacklisted_protein_coding$gene_id_base |
           Gene_name %in% blacklisted_protein_coding$gene_name)

cat(sprintf("  Blacklisted genes found in analysis: %d\n",
            nrow(blacklisted_in_analysis)))

if (nrow(blacklisted_in_analysis) > 0) {
  cat("\n  These genes will be flagged as 'High Genomic Repeats':\n")
  cat("   ", paste(blacklisted_in_analysis$Gene_name, collapse = ", "), "\n")
}

# =============================================================================
# STEP 6: Save results
# =============================================================================

cat("\nStep 6: Saving results...\n")

# Save full blacklist overlap results
blacklisted_genes %>%
  arrange(desc(overlap_fraction)) %>%
  write_csv("data/processed/blacklist_gene_overlaps.csv")
cat("  Saved: data/processed/blacklist_gene_overlaps.csv\n")

# Save list of blacklisted gene names for categorization script
# Include both protein-coding with significant overlap AND any genes in analysis
blacklisted_gene_names <- unique(c(
  blacklisted_protein_coding$gene_name,
  blacklisted_in_analysis$Gene_name
))

# If no genes found, also check genes with ANY overlap in the analysis set
if (length(blacklisted_gene_names) == 0) {
  cat("\n  No protein-coding genes with >= 10% overlap.\n")
  cat("  Checking for ANY overlap with analysis genes...\n")

  # Get all genes with any overlap
  all_blacklisted_names <- unique(blacklisted_genes$gene_name)

  # Match to analysis genes
  blacklisted_in_analysis <- chr21_analysis_genes %>%
    filter(gene_id_base %in% blacklisted_genes$gene_id_base |
             Gene_name %in% all_blacklisted_names)

  blacklisted_gene_names <- blacklisted_in_analysis$Gene_name

  cat(sprintf("  Analysis genes with any blacklist overlap: %d\n",
              length(blacklisted_gene_names)))
}

# Create simple lookup table (even if empty)
if (length(blacklisted_gene_names) > 0) {
  blacklist_lookup <- data.frame(
    Gene_name = blacklisted_gene_names,
    in_blacklist = TRUE,
    stringsAsFactors = FALSE
  )
} else {
  # Create empty dataframe with correct structure
  blacklist_lookup <- data.frame(
    Gene_name = character(0),
    in_blacklist = logical(0),
    stringsAsFactors = FALSE
  )
  cat("\n  NOTE: No chr21 protein-coding genes overlap blacklist regions.\n")
  cat("  This is expected - most blacklist regions are in centromeric areas\n")
  cat("  with few protein-coding genes.\n")
}

blacklist_lookup %>%
  write_csv("data/processed/blacklisted_genes.csv")
cat("  Saved: data/processed/blacklisted_genes.csv\n")

cat(sprintf("\n  Total blacklisted genes for categorization: %d\n",
            nrow(blacklist_lookup)))

# =============================================================================
# STEP 7: Summary
# =============================================================================

cat("\n=== Blacklist Processing Summary ===\n")
cat(sprintf("ENCODE blacklist regions on chr21: %d\n", nrow(blacklist_chr21)))
cat(sprintf("Total chr21 genes overlapping blacklist: %d\n",
            nrow(blacklisted_genes)))
cat(sprintf("Protein-coding genes with significant overlap: %d\n",
            nrow(blacklisted_protein_coding)))
cat(sprintf("Genes in current analysis to flag: %d\n",
            nrow(blacklist_lookup)))

if (nrow(blacklist_lookup) > 0) {
  cat("\nBlacklisted genes:\n")
  cat(" ", paste(blacklist_lookup$Gene_name, collapse = ", "), "\n")
} else {
  cat("\nNo protein-coding genes in the analysis overlap blacklist regions.\n")
  cat("The 'High Genomic Repeats' category may need alternative identification\n")
  cat("(e.g., genes with known paralogs or high repeat content from other sources).\n")
}

cat("\n=== Processing Complete ===\n")
cat("The file data/processed/blacklisted_genes.csv can be used by\n")
cat("02_categorize_genes.R to flag genes in blacklist regions.\n\n")
