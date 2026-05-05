# install_packages.R
#
# Purpose: Install all required R packages for T21-eQTL analysis
# Usage: Rscript install_packages.R
#
# This script will:
#   1. Check R version
#   2. Install CRAN packages
#   3. Install Bioconductor packages
#   4. Verify installations
#   5. Create renv.lock for reproducibility (optional)
#
# Author: Claude Code
# Date: 2025-11-11

cat("=== T21-eQTL Package Installation ===\n\n")

# =============================================================================
# STEP 1: Check R version
# =============================================================================

cat("Step 1: Checking R version...\n")
r_version <- R.version.string
cat(sprintf("  %s\n", r_version))

r_version_num <- as.numeric(R.version$major) +
                 as.numeric(R.version$minor) / 10

if (r_version_num < 4.0) {
  warning("R version < 4.0 detected. Some packages may not work properly.")
  cat("  Recommended: R >= 4.0\n")
} else {
  cat("  R version OK\n")
}

# =============================================================================
# STEP 2: Define required packages
# =============================================================================

cat("\nStep 2: Defining required packages...\n")

# CRAN packages
cran_packages <- c(
  # Data manipulation
  "tidyverse",      # Collection: dplyr, ggplot2, tidyr, readr, etc.
  "data.table",     # Fast data processing for large files

  # Visualization
  "ggplot2",        # Plotting (included in tidyverse but listed explicitly)
  "ggrepel",        # Label positioning for plots
  "ggalluvial",     # Alluvial/Sankey diagrams
  "RColorBrewer",   # Color palettes
  "pheatmap",       # Heatmaps
  "viridis",        # Perceptually uniform color scales

  # Utilities
  "here"            # Path management (optional but recommended)
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",         # Differential expression analysis
  "BiocManager"     # Bioconductor package manager
)

cat(sprintf("  CRAN packages: %d\n", length(cran_packages)))
cat(sprintf("  Bioconductor packages: %d\n", length(bioc_packages)))

# =============================================================================
# STEP 3: Install BiocManager (required for Bioconductor)
# =============================================================================

cat("\nStep 3: Installing BiocManager...\n")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  cat("  BiocManager installed\n")
} else {
  cat("  BiocManager already installed\n")
}

library(BiocManager)

# =============================================================================
# STEP 4: Install CRAN packages
# =============================================================================

cat("\nStep 4: Installing CRAN packages...\n")

for (pkg in cran_packages) {
  cat(sprintf("  Checking %s...", pkg))

  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(" installing...")
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org",
                      dependencies = TRUE)
      cat(" SUCCESS\n")
    }, error = function(e) {
      cat(" FAILED\n")
      cat(sprintf("    Error: %s\n", e$message))
    })
  } else {
    cat(" already installed\n")
  }
}

# =============================================================================
# STEP 5: Install Bioconductor packages
# =============================================================================

cat("\nStep 5: Installing Bioconductor packages...\n")

for (pkg in bioc_packages) {
  if (pkg == "BiocManager") next  # Already installed

  cat(sprintf("  Checking %s...", pkg))

  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(" installing...")
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      cat(" SUCCESS\n")
    }, error = function(e) {
      cat(" FAILED\n")
      cat(sprintf("    Error: %s\n", e$message))
    })
  } else {
    cat(" already installed\n")
  }
}

# =============================================================================
# STEP 6: Verify installations
# =============================================================================

cat("\nStep 6: Verifying installations...\n")

all_packages <- c(cran_packages, bioc_packages)
failed_packages <- character()

for (pkg in all_packages) {
  can_load <- require(pkg, character.only = TRUE, quietly = TRUE)

  if (!can_load) {
    failed_packages <- c(failed_packages, pkg)
    cat(sprintf("  FAILED: %s\n", pkg))
  }
}

if (length(failed_packages) == 0) {
  cat("  All packages installed successfully!\n")
} else {
  cat(sprintf("\n  WARNING: %d packages failed to install:\n",
              length(failed_packages)))
  for (pkg in failed_packages) {
    cat(sprintf("    - %s\n", pkg))
  }
  cat("\n  Try installing failed packages manually:\n")
  cat(sprintf("    install.packages(c(%s))\n",
              paste0("\"", paste(failed_packages, collapse = "\", \""), "\"")))
}

# =============================================================================
# STEP 7: Print package versions
# =============================================================================

cat("\nStep 7: Package versions:\n")

for (pkg in all_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(sprintf("  %-15s %s\n", pkg, version))
  }
}

# =============================================================================
# STEP 8: Optional - Initialize renv for reproducibility
# =============================================================================

cat("\nStep 8: Optional renv initialization...\n")
cat("  To create reproducible environment, run:\n")
cat("    install.packages('renv')\n")
cat("    renv::init()\n")
cat("    renv::snapshot()\n")
cat("\n  This creates renv.lock file tracking all package versions.\n")

# =============================================================================
# STEP 9: Save session info
# =============================================================================

cat("\nStep 9: Saving session info...\n")

if (!dir.exists("docs")) {
  dir.create("docs", recursive = TRUE)
}

writeLines(capture.output(sessionInfo()),
           "docs/package_installation_info.txt")
cat("  Saved: docs/package_installation_info.txt\n")

# =============================================================================
# Summary
# =============================================================================

cat("\n=== Installation Complete ===\n")

if (length(failed_packages) == 0) {
  cat("All packages installed successfully!\n")
  cat("You can now run the analysis scripts.\n")
} else {
  cat("Some packages failed to install.\n")
  cat("Please install them manually before running analysis.\n")
}

cat("\nNext steps:\n")
cat("  1. Run: Rscript scripts/00_preprocess_data.R\n")
cat("  2. Or run complete pipeline: bash scripts/run_all.sh\n\n")
