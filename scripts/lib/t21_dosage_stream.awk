# t21_dosage_stream.awk
#
# Stream a VCF-style PASS CSV (data/chr21_ds_PASS.csv) to a compact
# ALT-dosage matrix: one row per biallelic variant whose file-wide minor
# allele frequency (INFO AC/AN) is at least -v minmaf, one column per sample.
# Dosage is the count of non-reference alleles in the GT field (0-3 under
# trisomy); a GT containing "." is written empty (missing).
#
# Rows whose field count differs from the header's are skipped: a comma
# inside ALT or INFO marks a multi-allelic site, which the pipeline drops
# everywhere else too.
#
# Usage: awk -v minmaf=0.02 -f scripts/lib/t21_dosage_stream.awk data/chr21_ds_PASS.csv

BEGIN { FS = ","; OFS = "," }

NR == 1 {
  nf = NF
  printf "POS,REF,ALT"
  for (i = 10; i <= NF; i++) printf ",%s", $i
  printf "\n"
  next
}

NF != nf { next }

{
  ac = $8; an = $8
  if (!match(ac, /(^|;)AC=[0-9]+/)) next
  ac = substr(ac, RSTART, RLENGTH); sub(/.*AC=/, "", ac)
  if (!match(an, /(^|;)AN=[0-9]+/)) next
  an = substr(an, RSTART, RLENGTH); sub(/.*AN=/, "", an)
  if (an + 0 == 0) next
  af  = ac / an
  maf = (af < 0.5) ? af : 1 - af
  if (maf < minmaf) next

  line = $2 OFS $4 OFS $5
  for (i = 10; i <= NF; i++) {
    g = $i
    sub(/:.*/, "", g)
    if (index(g, ".")) line = line OFS ""
    else { n = gsub(/[1-9]/, "", g); line = line OFS n }
  }
  print line
}
