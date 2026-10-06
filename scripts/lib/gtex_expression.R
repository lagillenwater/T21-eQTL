# gtex_expression.R
#
# GTEx v10 gene expression across tissues, used to ask whether a chr21 gene is
# meaningfully expressed in whole blood relative to other tissues
# (manuscript issue #30). Same release as the eQTL data (v10).
#
# Source: the per-tissue median TPM table GTEx distributes with each release,
# one row per gene (GENCODE v39 ids, versioned) and one column per tissue
# (68 in v10). It is public (no requester-pays, no auth) and 8.8 MB, so it is
# downloaded once into data/gtex/ and read from there.
#   https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
#   gs://adult-gtex/bulk-gex/v10/rna-seq/

GTEX_MEDIAN_TPM_URL <- paste0(
  "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/",
  "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
GTEX_MEDIAN_TPM_FILE <- file.path("data", "gtex", basename(GTEX_MEDIAN_TPM_URL))
GTEX_WHOLE_BLOOD <- "Whole_Blood"

#' Read a GTEx median-TPM GCT file (plain or gzipped).
#'
#' Returns one row per gene: `ensembl_stable` (version stripped),
#' `gene_name`, then one numeric column per tissue. The `_PAR_Y` rows (the Y
#' copies of pseudoautosomal genes, which would repeat the X copy's stable id)
#' are dropped. Stops unless the file is GCT (first line `#1.2`) and carries a
#' Whole_Blood column, so a truncated or error-page download is never read as
#' data.
read_gtex_median_tpm <- function(path) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  lines <- readLines(con)
  if (!length(lines) || !startsWith(lines[1], "#1.2"))
    stop(sprintf("%s is not a GCT file (first line is not #1.2)", path), call. = FALSE)
  dt <- data.table::fread(text = lines[-(1:2)], sep = "\t", header = TRUE)
  if (!GTEX_WHOLE_BLOOD %in% names(dt))
    stop(sprintf("%s has no %s column", path, GTEX_WHOLE_BLOOD), call. = FALSE)
  dt <- dt[!endsWith(Name, "_PAR_Y")]
  dt[, Name := sub("[.].*$", "", Name)]
  data.table::setnames(dt, c("Name", "Description"), c("ensembl_stable", "gene_name"))
  dt[]
}

#' Path to the cached GTEx v10 median-TPM file, downloading it first if absent.
#'
#' The download goes to a `.part` file, is validated by reading it, and only
#' then takes the final name, so an interrupted fetch never leaves a file the
#' next run would trust.
gtex_median_tpm_file <- function(path = GTEX_MEDIAN_TPM_FILE, url = GTEX_MEDIAN_TPM_URL) {
  if (file.exists(path)) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  part <- paste0(path, ".part")
  cat(sprintf("  Downloading %s\n", url))
  ok <- tryCatch(utils::download.file(url, part, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) FALSE)
  if (!ok) stop(sprintf("could not download %s; fetch it manually to %s", url, path), call. = FALSE)
  read_gtex_median_tpm(part)
  file.rename(part, path)
  path
}

# Organ-system groups of the GTEx tissues, for display. Whole blood is its own
# group; the ileal lymphoid aggregate is lymphoid tissue, not gut.
GTEX_TISSUE_GROUPS <- c("Whole blood", "Immune and lymphoid", "Brain and nerve",
                        "Heart, vessels and muscle", "Digestive", "Reproductive", "Other")

#' Organ-system group of each GTEx tissue name (factor, GTEX_TISSUE_GROUPS order).
gtex_tissue_group <- function(tissue) {
  grp <- data.table::fcase(
    tissue == GTEX_WHOLE_BLOOD,                                      "Whole blood",
    grepl("^Spleen$|EBV-transformed_lymphocytes|Lymphode_Aggregate", tissue), "Immune and lymphoid",
    grepl("^Brain_|^Nerve_|^Pituitary$", tissue),                   "Brain and nerve",
    grepl("^Heart_|^Artery_|^Muscle_", tissue),                     "Heart, vessels and muscle",
    grepl("^Colon_|^Esophagus_|^Stomach|^Small_Intestine|^Liver|^Pancreas|^Minor_Salivary_Gland$", tissue),
                                                                    "Digestive",
    grepl("^Testis$|^Ovary$|^Uterus$|^Vagina$|^Cervix_|^Fallopian_Tube$|^Prostate$", tissue), "Reproductive",
    default = "Other")
  factor(grp, levels = GTEX_TISSUE_GROUPS)
}

#' Per-gene expression of one focal tissue against all tissues.
#'
#' `tpm`: one row per gene, `ensembl_stable` plus numeric tissue columns (any
#' other non-numeric column is ignored). Returns, per gene:
#'   focal_tpm          median TPM in the focal tissue
#'   focal_rank         its rank among the tissues, 1 = highest, ties averaged;
#'                      NA for a gene expressed in no tissue
#'   n_tissues          number of tissues ranked over
#'   top_tissue, top_tpm   the highest tissue and its TPM (top_tissue NA if none)
#'   focal_frac_of_top  focal TPM as a fraction of the top tissue's
#'   tau                tissue specificity (Yanai et al. 2005) on log2(TPM + 1):
#'                      0 = equal in every tissue, 1 = expressed in one only
tissue_expression_summary <- function(tpm, focal = GTEX_WHOLE_BLOOD) {
  tissues <- names(tpm)[vapply(tpm, is.numeric, logical(1))]
  if (!focal %in% tissues) stop(sprintf("no %s column among the tissues", focal), call. = FALSE)
  m <- as.matrix(tpm[, ..tissues])
  f <- m[, focal]
  top <- apply(m, 1, max)
  expressed <- top > 0
  above <- rowSums(m > f)
  tied  <- rowSums(m == f)
  x <- log2(m + 1)
  tau <- rowSums(1 - x / apply(x, 1, max)) / (length(tissues) - 1)
  data.table::data.table(
    ensembl_stable    = tpm$ensembl_stable,
    focal_tpm         = f,
    focal_rank        = ifelse(expressed, above + (tied + 1) / 2, NA_real_),
    n_tissues         = length(tissues),
    top_tissue        = ifelse(expressed, tissues[max.col(m, ties.method = "first")], NA_character_),
    top_tpm           = top,
    focal_frac_of_top = ifelse(expressed, f / top, NA_real_),
    tau               = ifelse(expressed, tau, NA_real_))
}

# ---- HTP expression on the TPM scale ------------------------------------------
# GTEx reports median TPM, so the cohort's counts are put on the same scale:
# reads per kilobase of gene, then scaled so each sample sums to one million.
# Gene length is the union of the gene's exons (GENCODE v44 basic, the GTF in
# data/), the same collapsed-model idea GTEx's RNA-SeQC uses; GTEx's own model
# is GENCODE v39 comprehensive, so lengths agree closely but not exactly.

GENCODE_GTF <- file.path("data", "gencode.v44.basic.annotation.gtf.gz")

#' Exon rows of a GTF (plain or gzipped): `ensembl_stable`, `start`, `end`.
#' The `_PAR_Y` copies are dropped, as in read_gtex_median_tpm().
read_gtf_exons <- function(path = GENCODE_GTF) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  lines <- readLines(con)
  lines <- lines[grepl("\texon\t", lines, fixed = TRUE)]
  dt <- data.table::fread(text = lines, sep = "\t", header = FALSE, select = c(4, 5, 9),
                          col.names = c("start", "end", "attributes"))
  dt[, gene_id := sub('^.*gene_id "([^"]+)".*$', "\\1", attributes)]
  dt <- dt[!endsWith(gene_id, "_PAR_Y")]
  dt[, .(ensembl_stable = sub("[.].*$", "", gene_id), start = as.integer(start), end = as.integer(end))]
}

#' Total length (bp) of the union of each gene's exons, coordinates 1-based
#' inclusive. An exon starting inside or right after the running end of the
#' gene's earlier exons extends the current block; any gap starts a new one.
exon_union_length <- function(exons) {
  ex <- data.table::copy(exons)[order(ensembl_stable, start, end)]
  # fill 0, not -Inf: shift() keeps the integer type and would turn -Inf into
  # NA; GTF coordinates start at 1, so 0 is below every exon.
  ex[, prev_end := data.table::shift(cummax(end), fill = 0L), by = ensembl_stable]
  ex[, block := cumsum(start > prev_end), by = ensembl_stable]
  blocks <- ex[, .(len = max(end) - min(start) + 1), by = .(ensembl_stable, block)]
  blocks[, .(length_bp = sum(len)), by = ensembl_stable]
}

#' Counts (genes x samples, gene row names) to TPM, given named gene lengths in
#' bp. Genes without a length are dropped before scaling.
counts_to_tpm <- function(counts, length_bp) {
  keep <- rownames(counts) %in% names(length_bp)
  rate <- counts[keep, , drop = FALSE] / (length_bp[rownames(counts)[keep]] / 1000)
  sweep(rate, 2, colSums(rate), "/") * 1e6
}
