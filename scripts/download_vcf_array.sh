#!/bin/bash
# Download the INCLUDE/Kids First VCF manifest as a SLURM job array on Alpine.
#
# Each array task downloads a stripe of the manifest -- rows where
# (row - 1) %% <array size> == <task index> - 1 -- so the files spread evenly
# across tasks. The token is read from ~/.cavatica_token by drs_download.R, so
# it is never written into this script or the job environment.
#
# Prerequisites (once, in the R that `module load` below provides):
#   Rscript -e 'install.packages(c("httr2","optparse","data.table"), repos="https://cloud.r-project.org")'
#   printf '%s' '<your Cavatica auth token>' > ~/.cavatica_token && chmod 600 ~/.cavatica_token
#
# Submit (edit the config block below first):
#   mkdir -p logs
#   sbatch scripts/download_vcf_array.sh /path/to/manifest_20260713_094446.csv
#
# The array size (--array=1-N) IS the number of chunks; the script reads it
# back from SLURM_ARRAY_TASK_COUNT, so to use more/fewer tasks just change the
# number after the dash. Keep the range 1-based and contiguous. The %K suffix
# caps how many tasks run at once (bandwidth / rate-limit throttle).

#SBATCH --job-name=vcf-drs
#SBATCH --partition=amilan
#SBATCH --qos=normal
##SBATCH --account=<your-allocation>   # uncomment/set if you have no default
#SBATCH --array=1-20%5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/vcf-drs-%A_%a.out
#SBATCH --error=logs/vcf-drs-%A_%a.err

set -euo pipefail

# ---- config: edit for your Alpine setup -----------------------------------
DRS_SCRIPT="${DRS_SCRIPT:-$HOME/T21-eQTL/scripts/drs_download.R}"
OUTDIR="${OUTDIR:-/pl/active/pivlab/projects/msubirana/dosage_comp/data/vcf}"
R_MODULE="${R_MODULE:-R}"     # R module name on Alpine (e.g. R/4.4.0)
# ---------------------------------------------------------------------------

MANIFEST="${1:?usage: sbatch scripts/download_vcf_array.sh <manifest.csv>}"

module purge
module load "$R_MODULE"

NCHUNKS="${SLURM_ARRAY_TASK_COUNT:?not a job array; submit with sbatch --array=1-N}"
CHUNK="${SLURM_ARRAY_TASK_ID}"

echo "task ${CHUNK}/${NCHUNKS} on $(hostname) -> ${OUTDIR}"

Rscript "$DRS_SCRIPT" "$MANIFEST" \
  --host    cavatica \
  --outdir  "$OUTDIR" \
  --workers "${SLURM_CPUS_PER_TASK:-4}" \
  --chunks  "$NCHUNKS" \
  --chunk   "$CHUNK"

echo "task ${CHUNK}/${NCHUNKS} done"
