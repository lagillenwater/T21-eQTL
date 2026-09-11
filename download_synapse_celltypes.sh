#!/usr/bin/env bash
# Download the HTP cell-type data from Synapse (folder syn31488783).
#
# Output: data/synapse/syn31488783/          every file in the Synapse folder,
#                                            mirrored with its own name
#         data/synapse/syn31488783_manifest.tsv  the Synapse manifest written by
#                                            the client (id, name, path, md5)
#         data/ is gitignored.
#
# Client. Uses the Synapse Python client (synapseclient, pinned below) in a
# virtual environment created under ~/.virtualenvs/synapser on first run. No
# R package is involved: the synapser R package does not currently install
# against this machine's R library (it pins rjson <= 0.2.21), and the Python
# client is what synapser wraps anyway.
#
# Authentication. A Synapse personal access token with the "download" scope
# is required (Synapse > Account Settings > Personal Access Tokens). Provide
# it one of two ways; it is never written into this script or the repo:
#   SYNAPSE_AUTH_TOKEN=<token>     environment variable, or
#   ~/.synapse_token               a file holding just the token on its first
#                                  line (chmod 600 it), or
#   ~/.synapseConfig               the Synapse client's own config file with an
#                                  [authentication] section (read by the client
#                                  itself; nothing is passed on the command line)
# Checked in that order; the first one found is used.
#
# Usage (from the repo root):
#   printf '%s' '<token>' > ~/.synapse_token && chmod 600 ~/.synapse_token
#   bash download_synapse_celltypes.sh
set -euo pipefail

SYN_ID="syn31488783"
DEST_DIR="data/synapse/${SYN_ID}"
MANIFEST="data/synapse/${SYN_ID}_manifest.tsv"
VENV="${SYNAPSE_VENV:-$HOME/.virtualenvs/synapser}"
CLIENT_VERSION="4.4.2"

# --- token ------------------------------------------------------------------
token="${SYNAPSE_AUTH_TOKEN:-}"
if [ -z "${token}" ] && [ -f "$HOME/.synapse_token" ]; then
  token="$(head -n 1 "$HOME/.synapse_token" | tr -d '[:space:]')"
fi
auth_args=()
if [ -n "${token}" ]; then
  auth_args=(-p "${token}")
elif [ -f "$HOME/.synapseConfig" ]; then
  echo "Using credentials from ~/.synapseConfig"
else
  cat >&2 <<MSG
No Synapse credentials found. Set SYNAPSE_AUTH_TOKEN, write the token to
~/.synapse_token (first line only, chmod 600), or place a .synapseConfig in
your home directory:
  printf '%s' '<token>' > ~/.synapse_token && chmod 600 ~/.synapse_token
MSG
  exit 1
fi

# --- client -----------------------------------------------------------------
if [ ! -x "${VENV}/bin/synapse" ]; then
  echo "Creating ${VENV} and installing synapseclient ${CLIENT_VERSION}..."
  python3 -m venv "${VENV}"
  "${VENV}/bin/pip" install --quiet --upgrade pip
  "${VENV}/bin/pip" install --quiet "synapseclient==${CLIENT_VERSION}"
fi
echo "synapse client: $("${VENV}/bin/synapse" --version 2>/dev/null | head -n 1)"

# --- download ---------------------------------------------------------------
mkdir -p "${DEST_DIR}"
echo "Syncing ${SYN_ID} -> ${DEST_DIR}"
# -r: recursive; --manifest: write the client's manifest beside the files;
# a token, when given, is passed to the client only (never exported).
"${VENV}/bin/synapse" "${auth_args[@]}" get -r "${SYN_ID}" \
  --downloadLocation "${DEST_DIR}" --manifest all

# The client names its manifest SYNAPSE_METADATA_MANIFEST.tsv inside the
# download location; keep a copy beside the folder under a stable name.
if [ -f "${DEST_DIR}/SYNAPSE_METADATA_MANIFEST.tsv" ]; then
  cp "${DEST_DIR}/SYNAPSE_METADATA_MANIFEST.tsv" "${MANIFEST}"
fi

n_files="$(find "${DEST_DIR}" -type f ! -name 'SYNAPSE_METADATA_MANIFEST.tsv' | wc -l | tr -d ' ')"
if [ "${n_files}" -eq 0 ]; then
  echo "No files were downloaded from ${SYN_ID}" >&2
  exit 1
fi

echo
echo "Downloaded ${n_files} file(s) to ${DEST_DIR}:"
find "${DEST_DIR}" -type f ! -name 'SYNAPSE_METADATA_MANIFEST.tsv' -exec ls -lh {} \; \
  | awk '{printf "  %8s  %s\n", $5, $NF}'
[ -f "${MANIFEST}" ] && echo "Manifest: ${MANIFEST}"
