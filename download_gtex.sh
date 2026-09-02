#!/usr/bin/env bash
# Download the GTEx v10 whole-blood cis-eQTL all-associations file for chr21.
#
# Output: data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet
# Consumer: scripts/02_filter_genotypes.R reads it with arrow::read_parquet()
#           and uses the columns gene_id, variant_id and pval_nominal.
#
# GTEx distributes the v10 all-associations results per tissue and per
# chromosome as parquet files in a requester-pays Google Cloud bucket
# (https://gtexportal.org/home/downloads/adult-gtex/qtl):
#   gs://gtex-resources/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/
# The file is used exactly as distributed; there is no conversion step, so no
# python or pyarrow is needed. Egress is billed to the Google Cloud project
# named in GCP_BILLING_PROJECT.
#
# Usage (from the repo root, with gcloud or gsutil installed and authenticated):
#   GCP_BILLING_PROJECT=<your-gcp-project> bash download_gtex.sh
set -euo pipefail

SRC_BUCKET="gs://gtex-resources"
SRC_DIR="GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations"
SRC_FILE="Whole_Blood.v10.allpairs.chr21.parquet"
SRC="${SRC_BUCKET}/${SRC_DIR}/${SRC_FILE}"
DEST="data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet"
MIN_BYTES=10000000   # the chr21 file is ~50 MB; anything far smaller is not the data

if [ -z "${GCP_BILLING_PROJECT:-}" ]; then
  cat >&2 <<MSG
GCP_BILLING_PROJECT is not set. The GTEx bucket is requester-pays, so a Google
Cloud project must be named to bill the egress:
  GCP_BILLING_PROJECT=<your-gcp-project> bash download_gtex.sh
MSG
  exit 1
fi

mkdir -p data
tmp="${DEST}.part"
rm -f "${tmp}"

if command -v gcloud >/dev/null 2>&1; then
  echo "Downloading ${SRC} with gcloud storage (billing project: ${GCP_BILLING_PROJECT})..."
  gcloud storage cp --billing-project="${GCP_BILLING_PROJECT}" "${SRC}" "${tmp}"
elif command -v gsutil >/dev/null 2>&1; then
  echo "Downloading ${SRC} with gsutil (billing project: ${GCP_BILLING_PROJECT})..."
  gsutil -u "${GCP_BILLING_PROJECT}" cp "${SRC}" "${tmp}"
else
  cat >&2 <<MSG
Neither gcloud nor gsutil is on PATH. Install the Google Cloud SDK
(https://cloud.google.com/sdk/docs/install), run 'gcloud auth login', and
re-run this script. Alternatively download the object by hand
  ${SRC}
(requester-pays: https://cloud.google.com/storage/docs/using-requester-pays)
and save it as
  ${DEST}
MSG
  exit 1
fi

size=$(stat -f%z "${tmp}" 2>/dev/null || stat -c%s "${tmp}")
if [ "${size}" -lt "${MIN_BYTES}" ]; then
  echo "Downloaded file is only ${size} bytes; expected tens of MB. Removing it." >&2
  rm -f "${tmp}"
  exit 1
fi
mv "${tmp}" "${DEST}"

# Confirm the file is what script 02 expects: a parquet carrying the columns it reads.
if command -v Rscript >/dev/null 2>&1; then
  Rscript --vanilla -e '
    suppressPackageStartupMessages(library(arrow))
    f <- commandArgs(trailingOnly = TRUE)[1]
    d <- read_parquet(f, as_data_frame = FALSE)
    need <- c("gene_id", "variant_id", "pval_nominal")
    miss <- setdiff(need, names(d))
    if (length(miss) > 0) stop("parquet is missing column(s): ", paste(miss, collapse = ", "))
    cat(sprintf("Parquet check OK: %d rows; columns: %s\n",
                d$num_rows, paste(names(d), collapse = ", ")))
  ' "${DEST}"
else
  echo "Rscript not found; skipping the parquet column check." >&2
fi

echo "Wrote ${DEST} (${size} bytes)"
