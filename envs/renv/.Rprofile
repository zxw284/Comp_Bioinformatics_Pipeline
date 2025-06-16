# R environment profile for bioinformatics project
# This file is sourced when R starts up in this directory

# Initialize renv if not already done
if (file.exists('renv/activate.R')) {
  source('renv/activate.R')
} else {
  # Instructions for first-time setup
  cat("\n=== R Environment Setup ===")
  cat("\nTo initialize renv for this project, run:")
  cat("\n  renv::init()")
  cat("\n  renv::install(c(")
  cat("\n    'BiocManager', 'Seurat', 'SingleCellExperiment',")
  cat("\n    'DESeq2', 'edgeR', 'limma', 'ComplexHeatmap',")
  cat("\n    'clusterProfiler', 'GSVA', 'fgsea',")
  cat("\n    'tidyverse', 'ggplot2', 'pheatmap'")
  cat("\n  ))")
  cat("\n  renv::snapshot()")
  cat("\n========================\n\n")
}

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Bioconductor setup
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
options(BioC_mirror = "https://bioconductor.org")

# Common project settings
options(
  stringsAsFactors = FALSE,
  max.print = 100,
  scipen = 10,
  width = 120
)

# Memory management for large datasets
if (!interactive()) {
  options(future.globals.maxSize = 8 * 1024^3)  # 8GB
}

