# gnomad.R
#
# ALT allele frequencies for chr21 positions from the gnomAD v4.1 genomes
# sites VCF, used as the independent check on which allele of a cis variant is
# the minor allele in the typical population (scripts/lib/alleles.R).
#
# The sites VCF is 7.8 GB, so it is never downloaded. It is public (no
# requester-pays, no auth) and tabix-indexed, and the index is 36 KB: a range
# query transfers only the compressed blocks covering the requested positions.
# Rsamtools::scanTabix does that over https.
#
# The result is cached under data/gnomad/ keyed on the positions queried, so
# the pipeline hits the network once per variant universe. Positions that
# gnomAD does not carry are recorded in the queried-positions file, so a gene
# whose variants are absent from gnomAD is not re-queried on every run.
#
# Provenance of the source file:
#   https://gnomad.broadinstitute.org/downloads#v4  (v4.1 genomes, GRCh38)
#   gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/
# GRCh38 throughout, matching GTEx v10's b38 variant ids.

GNOMAD_CHR21_URL <- paste0(
  "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/",
  "genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz")

# AF      global ALT allele frequency (all gnomAD genomes)
# AF_nfe  non-Finnish European, the closest single group to the GTEx whole
#         blood donors and to the HTP cohort; carried so a global-vs-European
#         disagreement about which allele is minor is visible rather than
#         hidden inside one number.
GNOMAD_AF_FIELDS <- c("AF", "AF_nfe")

GNOMAD_CACHE_DIR       <- file.path("data", "gnomad")
GNOMAD_CACHE_AF        <- file.path(GNOMAD_CACHE_DIR, "gnomad_chr21_af.csv")
GNOMAD_CACHE_POSITIONS <- file.path(GNOMAD_CACHE_DIR, "gnomad_chr21_queried_positions.txt")

#' Collapse positions into query ranges, merging any two within `gap` bp.
#'
#' One range per position costs one https round trip each (~0.3 s); merging
#' neighbours trades a few more returned lines for far fewer round trips.
#' @return data.table(start, end), sorted, non-overlapping.
gnomad_regions <- function(pos, gap = 2000L) {
  p <- sort(unique(as.numeric(pos[!is.na(pos)])))
  if (!length(p)) return(data.table::data.table(start = numeric(0), end = numeric(0)))
  grp <- cumsum(c(TRUE, diff(p) > gap))
  data.table::data.table(start = as.numeric(tapply(p, grp, min)),
                         end   = as.numeric(tapply(p, grp, max)))
}

#' Pull one INFO key out of VCF INFO strings, vectorized.
#' Anchored on a field boundary so AF does not match AF_nfe or AF_XX.
gnomad_info_field <- function(info, key) {
  m   <- regexpr(sprintf("(?<=^|;)%s=\\K[^;]*", key), info, perl = TRUE)
  out <- rep(NA_character_, length(info))
  ok  <- m > 0
  if (any(ok)) out[ok] <- regmatches(info, m)
  out
}

#' Parse sites-VCF data lines into one row per (POS, REF, ALT).
#'
#' gnomAD distributes one ALT per line, so no multi-allelic splitting is
#' needed; a line that still carries a comma in ALT is dropped rather than
#' guessed at. Header lines (starting '#') are ignored.
#' @return data.table(POS, REF, ALT, filter, one numeric column per field)
parse_gnomad_lines <- function(lines, fields = GNOMAD_AF_FIELDS) {
  keep <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(keep))
    return(cbind(data.table::data.table(POS = integer(0), REF = character(0),
                                        ALT = character(0), filter = character(0)),
                 stats::setNames(data.table::as.data.table(
                   rep(list(numeric(0)), length(fields))), fields)))
  f <- data.table::tstrsplit(keep, "\t", fixed = TRUE)
  out <- data.table::data.table(POS = as.integer(f[[2]]), REF = f[[4]],
                                ALT = f[[5]], filter = f[[7]])
  for (key in fields)
    data.table::set(out, j = key,
                    value = suppressWarnings(as.numeric(gnomad_info_field(f[[8]], key))))
  out <- out[!grepl(",", ALT, fixed = TRUE)]
  out[]
}

#' Fetch ALT allele frequencies for `pos` from the remote sites VCF.
#'
#' Ranges are queried in chunks and filtered down to `pos` as they arrive, so
#' the ~250 gnomAD variants per kb that a merged range returns never
#' accumulate in memory.
#' @return data.table(POS, REF, ALT, filter, AF, AF_nfe); zero rows on failure,
#'   with a "gnomad_status" attribute describing why.
fetch_gnomad_af <- function(pos, url = GNOMAD_CHR21_URL, fields = GNOMAD_AF_FIELDS,
                            gap = 2000L, chunk = 25L, verbose = TRUE,
                            index_path = file.path(GNOMAD_CACHE_DIR,
                                                   basename(paste0(url, ".tbi")))) {
  empty <- parse_gnomad_lines(character(0), fields)
  if (!requireNamespace("Rsamtools", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    attr(empty, "gnomad_status") <- "Rsamtools/GenomicRanges not installed"
    return(empty)
  }
  regions <- gnomad_regions(pos, gap)
  if (!nrow(regions)) {
    attr(empty, "gnomad_status") <- "no positions requested"
    return(empty)
  }
  want  <- unique(as.integer(pos[!is.na(pos)]))
  idx   <- split(seq_len(nrow(regions)), ceiling(seq_len(nrow(regions)) / chunk))
  parts <- vector("list", length(idx))
  # Keep the index beside the cache. Left to itself htslib drops the .tbi in
  # the working directory, which for this pipeline is the repository root.
  dir.create(dirname(index_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(index_path)) {
    got <- try(utils::download.file(paste0(url, ".tbi"), index_path,
                                    mode = "wb", quiet = !verbose), silent = TRUE)
    if (inherits(got, "try-error")) {
      unlink(index_path)
      attr(empty, "gnomad_status") <- paste0("could not download the tabix index: ",
                                             trimws(as.character(got)))
      return(empty)
    }
  }
  tf    <- Rsamtools::TabixFile(url, index = index_path)
  t0    <- Sys.time()
  for (i in seq_along(idx)) {
    rows <- regions[idx[[i]]]
    gr   <- GenomicRanges::GRanges("chr21", IRanges::IRanges(start = rows$start, end = rows$end))
    res  <- try(Rsamtools::scanTabix(tf, param = gr), silent = TRUE)
    if (inherits(res, "try-error")) {
      attr(empty, "gnomad_status") <- paste0("tabix query failed: ",
                                             trimws(as.character(res)))
      return(empty)
    }
    parsed    <- parse_gnomad_lines(unlist(res, use.names = FALSE), fields)
    parts[[i]] <- parsed[POS %in% want]
    if (verbose && (i %% 4 == 0 || i == length(idx)))
      cat(sprintf("\r    gnomAD: %d / %d ranges (%.0f s)",
                  min(i * chunk, nrow(regions)), nrow(regions),
                  as.numeric(Sys.time() - t0, units = "secs")))
  }
  if (verbose) cat("\n")
  out <- unique(data.table::rbindlist(parts), by = c("POS", "REF", "ALT"))
  attr(out, "gnomad_status") <- "ok"
  out[]
}

#' Read the cached gnomAD frequencies and the positions already queried.
#' @return list(af = data.table, queried = integer vector)
read_gnomad_cache <- function(af_path = GNOMAD_CACHE_AF,
                              pos_path = GNOMAD_CACHE_POSITIONS,
                              fields = GNOMAD_AF_FIELDS) {
  af <- if (file.exists(af_path)) data.table::fread(af_path) else parse_gnomad_lines(character(0), fields)
  queried <- if (file.exists(pos_path)) as.integer(readLines(pos_path)) else integer(0)
  list(af = af, queried = queried[!is.na(queried)])
}

#' Cached gnomAD ALT frequencies for `pos`, fetching only what is missing.
#'
#' Offline or with the fetch disabled (T21_GNOMAD_OFFLINE=1) it returns
#' whatever the cache holds; the caller sees which positions are covered and
#' reports the gap rather than treating a missing frequency as a known one.
#' @return data.table with a "gnomad_status" attribute.
ensure_gnomad_af <- function(pos, af_path = GNOMAD_CACHE_AF,
                             pos_path = GNOMAD_CACHE_POSITIONS,
                             fields = GNOMAD_AF_FIELDS, verbose = TRUE, ...) {
  cache   <- read_gnomad_cache(af_path, pos_path, fields)
  want    <- unique(as.integer(pos[!is.na(pos)]))
  missing <- setdiff(want, cache$queried)
  offline <- nzchar(Sys.getenv("T21_GNOMAD_OFFLINE", ""))
  if (length(missing) && !offline) {
    if (verbose)
      cat(sprintf("  gnomAD: %d of %d positions not cached; querying %s\n",
                  length(missing), length(want), basename(GNOMAD_CHR21_URL)))
    fetched <- fetch_gnomad_af(missing, fields = fields, verbose = verbose, ...)
    if (identical(attr(fetched, "gnomad_status"), "ok")) {
      dir.create(dirname(af_path), recursive = TRUE, showWarnings = FALSE)
      cache$af      <- unique(data.table::rbindlist(list(cache$af, fetched), use.names = TRUE),
                              by = c("POS", "REF", "ALT"))
      cache$queried <- sort(union(cache$queried, missing))
      data.table::fwrite(cache$af, af_path)
      writeLines(as.character(cache$queried), pos_path)
    } else {
      out <- cache$af[POS %in% want]
      attr(out, "gnomad_status") <- attr(fetched, "gnomad_status")
      return(out[])
    }
  }
  out <- cache$af[POS %in% want]
  covered <- length(intersect(want, cache$queried))
  attr(out, "gnomad_status") <- if (covered == length(want)) "ok" else
    sprintf("%d of %d positions not queried%s", length(want) - covered, length(want),
            if (offline) " (T21_GNOMAD_OFFLINE set)" else "")
  out[]
}
