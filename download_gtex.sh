#!/bin/bash
# Script to manually download GTEx Whole Blood eQTL data

echo "Attempting to download GTEx v8 Whole Blood eQTL data..."

# Try the new GTEx URL structure
curl -L -o data/gtex_whole_blood_eqtls.txt.gz \
  "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"

if [ $? -eq 0 ] && [ -f data/gtex_whole_blood_eqtls.txt.gz ]; then
    filesize=$(stat -f%z data/gtex_whole_blood_eqtls.txt.gz 2>/dev/null || stat -c%s data/gtex_whole_blood_eqtls.txt.gz 2>/dev/null)
    if [ "$filesize" -gt 1000000 ]; then
        echo "Download successful! File size: $filesize bytes"
        exit 0
    else
        echo "Download failed (file too small)"
        rm -f data/gtex_whole_blood_eqtls.txt.gz
    fi
fi

echo "Trying alternative method with wget..."
wget -O data/gtex_whole_blood_eqtls.txt.gz \
  "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"

if [ $? -eq 0 ] && [ -f data/gtex_whole_blood_eqtls.txt.gz ]; then
    echo "Download successful with wget!"
    exit 0
fi

echo ""
echo "Automatic download failed. Please download manually:"
echo "1. Visit: https://gtexportal.org/home/datasets"
echo "2. Click on 'Download' -> 'Single-Tissue cis-QTL Data'"
echo "3. Find and download: Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"
echo "4. Save to: data/gtex_whole_blood_eqtls.txt.gz"
exit 1
