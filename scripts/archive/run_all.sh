#!/bin/bash
# run_all.sh
#
# Purpose: Run complete T21-eQTL analysis pipeline
# Usage: bash scripts/run_all.sh

echo "==================================="
echo "T21-eQTL Analysis Pipeline"
echo "==================================="
echo ""

# Check if we're in the right directory
if [ ! -f "scripts/00_preprocess_data.R" ]; then
    echo "Error: Must run from repository root directory"
    echo "Usage: bash scripts/run_all.sh"
    exit 1
fi

# Function to check if script succeeded
check_success() {
    if [ $? -eq 0 ]; then
        echo "SUCCESS: $1"
        echo ""
    else
        echo "ERROR: $1 failed"
        echo "Check error messages above"
        exit 1
    fi
}

# Step 0: Preprocess data
echo "Step 0: Preprocessing data (convert long to wide format)..."
Rscript scripts/00_preprocess_data.R
check_success "Data preprocessing"

# Step 1: DESeq2 analysis
echo "Step 1: Running trisomy-aware DESeq2 analysis..."
Rscript scripts/01_deseq2_analysis.R
check_success "DESeq2 analysis"

# Step 2: Categorize genes
echo "Step 2: Categorizing chr21 genes..."
Rscript scripts/02_categorize_genes.R
check_success "Gene categorization"

# Step 3: Volcano plot
echo "Step 3: Creating volcano plots..."
Rscript scripts/03_volcano_plot.R
check_success "Volcano plot"

# Step 4: Alluvial plot
echo "Step 4: Creating alluvial diagram (Panel D)..."
Rscript scripts/04_alluvial_plot.R
check_success "Alluvial plot"

# Step 5: eQTL analysis
echo "Step 5: Running eQTL cross-reference..."
echo "NOTE: This requires GTEx eQTL data (see script for details)"
Rscript scripts/05_eqtl_analysis.R
check_success "eQTL analysis"

echo "==================================="
echo "Pipeline Complete!"
echo "==================================="
echo ""
echo "Results saved to:"
echo "  - results/tables/"
echo "  - results/figures/"
echo ""
echo "Main outputs:"
echo "  - results/figures/panel_D_alluvial.pdf"
echo "  - results/tables/chr21_genes_categorized.csv"
echo "  - results/tables/final_gene_classification.csv"
echo ""
