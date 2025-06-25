#!/usr/bin/env Rscript
# Package Management Module
# Description: Handles installation and loading of required packages

# Function to install CRAN packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Function to install Bioconductor packages
install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Function to install GitHub packages
install_if_missing_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github(repo, dependencies = TRUE)
  }
}

# Define package sets
.get_package_sets <- function() {
  list(
    bioc_pkgs = c(
      "AnnotationHub", "biomaRt", "clusterProfiler", "decoupleR",
      "dorothea", "enrichplot", "ensembldb", "KEGGREST", "MAST",
      "OmnipathR", "scater", "viper", "zinbwave", "sparseMatrixStats",
      "progeny", "singscore", "GSVA", "GSEABase"
    ),
    cran_pkgs = c(
      "bigmemory", "cowplot", "data.table", "dplyr", "DT", "flexmix",
      "flextable", "foreach", "future", "ggplot2", "ggpubr", "igraph",
      "KernSmooth", "kableExtra", "knitr", "magrittr", "msigdbr",
      "patchwork", "pheatmap", "RCurl", "RColorBrewer", "ROCR", "SAVER",
      "sctransform", "Seurat", "stringi", "stringr", "tibble", "tidyr",
      "uwot", "writexl", "xlsx", "checkmate", "clustermole"
    ),
    gh_pkgs = c(
      "SeuratWrappers" = "satijalab/seurat-wrappers",
      "Stringendo" = "vertesy/Stringendo",
      "CodeAndRoll2" = "vertesy/CodeAndRoll2",
      "ReadWriter" = "vertesy/ReadWriter",
      "MarkdownHelpers" = "vertesy/MarkdownHelpers",
      "MarkdownReports" = "vertesy/MarkdownReports",
      "ggExpress" = "vertesy/ggExpress",
      "DatabaseLinke.R" = "vertesy/DatabaseLinke.R",
      "Seurat.utils" = "vertesy/Seurat.utils",
      "SeuratData" = "satijalab/seurat-data",
      "DoubletFinder" = "chris-mcginnis-ucsf/DoubletFinder",
      "Lamian" = "Winnie09/Lamian",
      "Libra" = "cran/Libra",
      "PISCES" = "califano-lab/PISCES",
      "SCPA" = "jackbibby1/SCPA",
      "colorblindr" = "clauswilke/colorblindr",
      "multicross" = "cran/multicross"
    )
  )
}

# Install all required packages
install_required_packages <- function() {
  # Set repository
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  options(repos = BiocManager::repositories())
  
  # Install r-spatial/s2
  remotes::install_github("r-spatial/s2@master", dependencies = TRUE)
  
  # Get package sets
  pkgs <- .get_package_sets()
  
  # Install packages
  lapply(pkgs$bioc_pkgs, install_if_missing_bioc)
  lapply(pkgs$cran_pkgs, install_if_missing)
  mapply(install_if_missing_github, 
         names(pkgs$gh_pkgs), 
         pkgs$gh_pkgs, 
         USE.NAMES = FALSE)
  
  # Install specific version of qpcR
  remotes::install_version("qpcR", version = "1.4-1",
                          repos = "https://cran.r-project.org")
}

# Load all required packages
load_required_packages <- function() {
  pkgs <- .get_package_sets()
  invisible(lapply(c(pkgs$bioc_pkgs, pkgs$cran_pkgs, names(pkgs$gh_pkgs)),
                  function(pkg) {
                    suppressPackageStartupMessages(
                      library(pkg, character.only = TRUE, quietly = TRUE)
                    )
                  }))
}

# Configure parallel processing if flexiblas is available
configure_parallel_processing <- function() {
  if (require("flexiblas", quietly = TRUE)) {
    flexiblas_load_backend("OPENBLAS-SERIAL")
    flexiblas_switch(n = 2)
    flexiblas_set_num_threads(1)
    message("Configured FlexiBLAS for parallel operations")
  } else {
    message("FlexiBLAS not available. For best performance, consider installing it.")
  }
}

# Main function to set up environment
setup_environment <- function() {
  message("Setting up environment...")
  install_required_packages()
  load_required_packages()
  configure_parallel_processing()
  message("Environment setup complete.")
}

# Run setup if script is executed directly
if (sys.nframe() == 0) {
  setup_environment()
}
