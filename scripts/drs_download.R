#!/usr/bin/env Rscript

# Download files listed in an INCLUDE / Kids First DRS manifest.
# R port of drs_download.py.
#
# Two ways to address each file:
#   cavatica  -> drs://cavatica-ga4gh-api.sbgenomics.com/<file id>   (Seven Bridges)
#               the file id is the manifest's `id` column (or a `drs_uri` column
#               if present); auth is the X-SBG-Auth-Token header.
#   gen3      -> the manifest's `access_url` column
#               (drs://nci-crdc.datacommons.io/dg.4DFC/...); auth is a Bearer token.
#
# Each host requires its own token. Provide it either as an environment
# variable or in a token file (env var wins if both are present):
#
#   CAVATICA_TOKEN  or  ~/.cavatica_token   Seven Bridges auth token
#                                           (Cavatica > Developer > Auth token)
#   GEN3_TOKEN      or  ~/.gen3_token       Gen3/DCF access token
#
# The token file should hold just the token on the first line; chmod 600 it.
# Override the file path with --token-file.
#
# Resolution follows the GA4GH DRS 1.x spec:
#   1. GET  {https_base}/ga4gh/drs/v1/objects/{object_id}
#   2. GET  .../objects/{object_id}/access/{access_id}   -> signed https URL
#   3. Stream the signed URL to disk.
#
# Usage:
#   printf '%s' '<token>' > ~/.cavatica_token && chmod 600 ~/.cavatica_token
#   Rscript scripts/drs_download.R manifest.csv --host cavatica
#   Rscript scripts/drs_download.R manifest.csv --host gen3 --workers 4
#   # one chunk of a job array (see scripts/download_vcf_array.sh):
#   Rscript scripts/drs_download.R manifest.csv --chunks 20 --chunk 3
#
# Signed URLs are short-lived, so resolution happens immediately before each
# download rather than in a separate up-front pass.
#
# Requires: httr2 (>= 1.0.0), optparse, data.table. Install with
#   install.packages(c("httr2", "optparse", "data.table"))

suppressPackageStartupMessages({
  have <- vapply(c("httr2", "optparse", "data.table"),
                 requireNamespace, logical(1), quietly = TRUE)
})
if (!all(have)) {
  stop("missing packages: install.packages(c('httr2','optparse','data.table'))",
       call. = FALSE)
}
library(httr2)
library(optparse)

CHUNK_KB <- 8L * 1024L   # 8 MiB streaming buffer
RETRIES  <- 4L
BACKOFF  <- 3            # seconds, doubled each retry
CAVATICA_DRS_HOST <- "cavatica-ga4gh-api.sbgenomics.com"
DEFAULT_OUTDIR <- "/pl/active/pivlab/projects/msubirana/dosage_comp/data/vcf"

## ---- helpers ---------------------------------------------------------------

# Read a manifest cell, coalescing NULL/NA/absent columns to "".
get_field <- function(row, key) {
  v <- row[[key]]
  if (is.null(v) || length(v) == 0 || is.na(v)) "" else trimws(as.character(v)[1])
}

# Extension of a filename, preserving compound ".xxx.gz" doubles (e.g. ".vcf.gz")
# and otherwise the final ".ext". Returns "" when there is no dot.
file_ext <- function(fn) {
  fn <- basename(fn)
  if (grepl("\\.[A-Za-z0-9]+\\.gz$", fn)) {
    return(sub(".*?(\\.[A-Za-z0-9]+\\.gz)$", "\\1", fn))
  }
  m <- regmatches(fn, regexpr("\\.[A-Za-z0-9]+$", fn))
  if (length(m)) m else ""
}

# Output filename for a row: <external_participant_id> + the source extension,
# so files are named by participant while staying tool-readable (.vcf.gz).
# Falls back to the source file name when no participant id is present.
out_filename <- function(row, src_name) {
  pid <- get_field(row, "external_participant_id")
  if (pid != "") paste0(pid, file_ext(src_name)) else src_name
}

# Split a drs:// URI into (https_base, object_id).
parse_drs <- function(uri) {
  if (is.null(uri) || !startsWith(uri, "drs://")) {
    stop(sprintf("not a DRS URI: %s", uri))
  }
  rest  <- sub("^drs://", "", uri)
  slash <- regexpr("/", rest, fixed = TRUE)
  if (slash < 1L) stop(sprintf("no object id in DRS URI: %s", uri))
  list(
    base      = paste0("https://", substr(rest, 1L, slash - 1L)),
    object_id = substr(rest, slash + 1L, nchar(rest))   # may itself contain '/'
  )
}

# Resolve a token: the environment variable wins; otherwise fall back to the
# first line of token_file (trimmed) if that file exists.
read_token <- function(env_name, token_file) {
  v <- Sys.getenv(env_name, unset = "")
  if (v != "") return(v)
  path <- path.expand(token_file)
  if (file.exists(path)) return(trimws(readLines(path, n = 1, warn = FALSE)))
  ""
}

# Auth header list keyed off the DRS host. Seven Bridges wants the token in an
# X-SBG-Auth-Token header; Gen3/DCF wants an Authorization Bearer token.
auth_header <- function(host, cavatica_token, gen3_token) {
  if (grepl("sbgenomics.com", host, fixed = TRUE)) {
    if (is.null(cavatica_token) || cavatica_token == "") {
      stop("CAVATICA_TOKEN is not set")
    }
    return(list("X-SBG-Auth-Token" = cavatica_token))
  }
  if (is.null(gen3_token) || gen3_token == "") stop("GEN3_TOKEN is not set")
  list(Authorization = paste("Bearer", gen3_token))
}

# Build the DRS URI for a manifest row given the chosen host. Cavatica prefers
# an explicit drs_uri column but falls back to constructing one from the file
# id; Gen3 uses the access_url column verbatim.
drs_uri_for_row <- function(row, host) {
  if (host == "gen3") return(get_field(row, "access_url"))
  u <- get_field(row, "drs_uri")
  if (u != "") return(u)
  fid <- get_field(row, "id")
  if (fid != "") return(paste0("drs://", CAVATICA_DRS_HOST, "/", fid))
  ""
}

# Resolve a DRS URI to a signed https URL plus a suggested name and size.
resolve <- function(drs_uri, auth, access_id = NULL) {
  parsed  <- parse_drs(drs_uri)
  obj_url <- paste0(parsed$base, "/ga4gh/drs/v1/objects/", parsed$object_id)

  obj <- request(obj_url) |>
    req_headers(!!!auth) |>
    req_timeout(60) |>
    req_perform() |>
    resp_body_json(check_type = FALSE)

  methods <- obj$access_methods
  if (is.null(methods) || length(methods) == 0) {
    stop(sprintf("no access_methods for %s", drs_uri))
  }

  has_url <- function(m) !is.null(m$access_url) && !is.null(m$access_url$url)

  # Prefer the requested access_id; otherwise take the first method that
  # either carries a direct access_url or an access_id we can exchange.
  chosen <- NULL
  if (!is.null(access_id)) {
    chosen <- Find(function(m) identical(m$access_id, access_id), methods)
  }
  if (is.null(chosen)) {
    chosen <- Find(function(m) has_url(m) || !is.null(m$access_id), methods)
  }
  if (is.null(chosen)) stop(sprintf("no usable access method for %s", drs_uri))

  if (has_url(chosen)) {
    signed <- chosen$access_url$url
  } else {
    aid <- chosen$access_id
    signed <- request(paste0(obj_url, "/access/", aid)) |>
      req_headers(!!!auth) |>
      req_timeout(60) |>
      req_perform() |>
      resp_body_json(check_type = FALSE) |>
      (\(x) x$url)()
  }

  list(url = signed, name = obj$name, size = obj$size)
}

# Stream a signed URL to disk. Resumes if a partial .part file exists.
download_file <- function(url, dest, expected_size = NULL) {
  part <- paste0(dest, ".part")
  have <- if (file.exists(part)) file.info(part)$size else 0

  # Signed URLs already embed credentials; do not send Authorization.
  # connecttimeout caps the handshake; low_speed_* aborts a stalled transfer
  # without capping total time (large genomics files may take a while).
  req <- request(url) |>
    req_options(connecttimeout = 30, low_speed_limit = 1, low_speed_time = 300) |>
    req_error(is_error = function(resp) {
      s <- resp_status(resp)
      s >= 400 && s != 416          # 416 == range already satisfied -> complete
    })

  mode <- "wb"
  if (have > 0) {
    req  <- req |> req_headers(Range = sprintf("bytes=%.0f-", have))
    mode <- "ab"
  }

  resp <- req_perform_connection(req)
  on.exit(close(resp), add = TRUE)

  if (resp_status(resp) != 416) {
    con_out <- file(part, open = mode)
    on.exit(close(con_out), add = TRUE)
    repeat {
      chunk <- resp_stream_raw(resp, kb = CHUNK_KB)
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }
  }

  actual <- file.info(part)$size
  if (!is.null(expected_size) && !is.na(expected_size) && actual != expected_size) {
    stop(sprintf("size mismatch for %s: got %.0f, expected %.0f",
                 basename(dest), actual, expected_size))
  }

  if (file.exists(dest)) unlink(dest)
  file.rename(part, dest)
  file.info(dest)$size
}

# Resolve + download one manifest row, with retry/backoff. Never throws;
# returns list(fname, status, detail).
handle_row <- function(row, host, outdir, cavatica_token, gen3_token) {
  drs_uri <- drs_uri_for_row(row, host)
  name    <- get_field(row, "file_name")
  if (name == "") name <- get_field(row, "name")

  if (drs_uri == "") {
    return(list(fname  = if (name != "") name else "<unnamed>",
                status = "skipped",
                detail = if (host == "cavatica") "no drs_uri or id in row"
                         else "no access_url in row"))
  }

  last_err <- NULL
  for (attempt in seq_len(RETRIES)) {
    result <- tryCatch({
      auth <- auth_header(parse_drs(drs_uri)$base, cavatica_token, gen3_token)
      r    <- resolve(drs_uri, auth)

      src_name <- if (name != "") {
        name
      } else if (!is.null(r$name) && nzchar(r$name)) {
        r$name
      } else {
        basename(url_parse(r$url)$path)
      }
      fname <- out_filename(row, src_name)
      dest  <- file.path(outdir, fname)

      if (file.exists(dest) && (is.null(r$size) || file.info(dest)$size == r$size)) {
        list(fname = fname, status = "exists", detail = "")
      } else {
        n <- download_file(r$url, dest, r$size)
        list(fname = fname, status = "ok", detail = sprintf("%.0f bytes", n))
      }
    }, error = function(e) {
      last_err <<- conditionMessage(e)
      NULL
    })

    if (!is.null(result)) return(result)
    if (attempt < RETRIES) Sys.sleep(BACKOFF * 2^(attempt - 1))
  }

  list(fname  = if (name != "") name else drs_uri,
       status = "failed",
       detail = if (is.null(last_err)) "unknown error" else last_err)
}

## ---- main ------------------------------------------------------------------

main <- function() {
  option_list <- list(
    make_option("--outdir",  type = "character", default = DEFAULT_OUTDIR,
                help = "output directory [default %default]"),
    make_option("--host",    type = "character", default = "cavatica",
                help = "cavatica -> resolve via the file id (or drs_uri) column; gen3 -> use the access_url column [default %default]"),
    make_option("--workers", type = "integer",   default = 4L,
                help = "parallel download workers [default %default]"),
    make_option("--token-file", type = "character", default = NA_character_,
                dest = "token_file",
                help = "path to a file holding the token (overrides ~/.cavatica_token / ~/.gen3_token)"),
    make_option("--chunks",  type = "integer",   default = NA_integer_,
                help = "split the manifest into this many chunks (e.g. the SLURM array size)"),
    make_option("--chunk",   type = "integer",   default = NA_integer_,
                help = "1-based chunk index to download (e.g. SLURM_ARRAY_TASK_ID)"),
    make_option("--limit",   type = "integer",   default = NA_integer_,
                help = "download only the first N rows (applied before chunking)")
  )
  parser <- OptionParser(usage = "%prog [options] manifest.csv",
                         option_list = option_list,
                         description = "Download files from a DRS manifest.")
  parsed   <- parse_args(parser, positional_arguments = 1)
  opt      <- parsed$options
  manifest <- parsed$args[1]

  if (!opt$host %in% c("cavatica", "gen3")) {
    stop("--host must be 'cavatica' or 'gen3'", call. = FALSE)
  }

  cav_file  <- if (!is.na(opt$token_file)) opt$token_file else "~/.cavatica_token"
  gen3_file <- if (!is.na(opt$token_file)) opt$token_file else "~/.gen3_token"
  cavatica_token <- read_token("CAVATICA_TOKEN", cav_file)
  gen3_token     <- read_token("GEN3_TOKEN",     gen3_file)
  if (opt$host == "cavatica" && cavatica_token == "") {
    stop(sprintf("no Cavatica token: set CAVATICA_TOKEN or put it in %s",
                 cav_file), call. = FALSE)
  }
  if (opt$host == "gen3" && gen3_token == "") {
    stop(sprintf("no Gen3 token: set GEN3_TOKEN or put it in %s",
                 gen3_file), call. = FALSE)
  }

  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  df <- data.table::fread(manifest, colClasses = "character",
                          na.strings = c("", "NA"), showProgress = FALSE)
  if (opt$host == "gen3") {
    if (!"access_url" %in% names(df)) {
      stop("column 'access_url' not found in manifest", call. = FALSE)
    }
    uri_src <- "access_url"
  } else {
    if (!any(c("drs_uri", "id") %in% names(df))) {
      stop("cavatica mode needs a 'drs_uri' or 'id' column", call. = FALSE)
    }
    uri_src <- if ("drs_uri" %in% names(df)) "drs_uri" else "id"
  }

  rows <- lapply(seq_len(nrow(df)), function(i) as.list(df[i]))
  if (!is.na(opt$limit)) rows <- utils::head(rows, opt$limit)

  # Chunk selection for job arrays: keep every (chunks)-th row starting at chunk.
  # Striping (modulo) rather than contiguous slices balances the load evenly.
  if (!is.na(opt$chunks)) {
    if (opt$chunks < 1L) stop("--chunks must be >= 1", call. = FALSE)
    if (is.na(opt$chunk) || opt$chunk < 1L || opt$chunk > opt$chunks) {
      stop("--chunk must be an integer in 1..chunks", call. = FALSE)
    }
    keep <- ((seq_along(rows) - 1L) %% opt$chunks) == (opt$chunk - 1L)
    rows <- rows[keep]
    message(sprintf("chunk %d of %d", opt$chunk, opt$chunks))
  }
  n <- length(rows)

  workers <- opt$workers
  if (.Platform$OS.type == "windows" && workers > 1L) {
    message("parallel downloads unsupported on Windows; using 1 worker")
    workers <- 1L
  }

  message(sprintf("%d rows; host=%s; uri source=%s; workers=%d",
                  n, opt$host, uri_src, workers))

  results <- parallel::mclapply(
    rows,
    function(r) handle_row(r, opt$host, opt$outdir, cavatica_token, gen3_token),
    mc.cores       = max(1L, workers),
    mc.preschedule = FALSE
  )

  counts <- c(ok = 0L, exists = 0L, failed = 0L, skipped = 0L)
  for (i in seq_len(n)) {
    res <- results[[i]]
    if (inherits(res, "try-error")) {   # worker died outside handle_row's tryCatch
      res <- list(fname = "<error>", status = "failed", detail = as.character(res))
    }
    counts[res$status] <- counts[res$status] + 1L
    message(sprintf("[%d/%d] %-7s %s %s", i, n, res$status, res$fname, res$detail))
  }

  message(sprintf("\ndone: %d downloaded, %d already present, %d failed, %d skipped",
                  counts[["ok"]], counts[["exists"]], counts[["failed"]], counts[["skipped"]]))
  if (counts[["failed"]] > 0L) quit(status = 1L, save = "no")
}

main()
