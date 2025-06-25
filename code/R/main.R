#!/usr/bin/env Rscript
# Main Analysis Script
# Description: Orchestrates the scRNA-seq analysis pipeline

# Set up logging
log_file <- file("logs/analysis_log.txt", open = "wt")
sink(file = log_file, type = "message")

# Load utility functions
source("src/load_packages.R")
source("src/load_data.R")
source("src/quality_control.R")
# Add other module sources as they are created
# source("src/normalization.R")
# source("src/dimensionality_reduction.R")
# source("src/clustering.R")
# source("src/differential_expression.R")
# source("src/visualization.R")
# source("src/trajectory_analysis.R")
# source("src/pathway_analysis.R")

# Main function to run the entire pipeline
run_analysis <- function(config) {
  # Create necessary directories
  dir.create("results", showWarnings = FALSE, recursive = TRUE)
  dir.create("logs", showWarnings = FALSE, recursive = TRUE)
  
  # Start timing
  start_time <- Sys.time()
  message("Starting analysis at: ", start_time)
  
  # 1. Set up environment
  message("\n=== Setting up environment ===")
  setup_environment()
  
  # 2. Load data
  message("\n=== Loading data ===")
  data <- load_scrnaseq_data(config)
  
  # 3. Quality control
  message("\n=== Performing quality control ===")
  seurat_obj <- perform_quality_control(
    data$seurat,
    cell_cycle_genes = data$cell_cycle_genes,
    species = config$species,
    min_genes = config$qc$min_genes,
    max_genes = config$qc$max_genes,
    max_mt = config$qc$max_mt,
    min_counts = config$qc$min_counts,
    max_counts = config$qc$max_counts
  )
  
  # Save session info
  capture.output(sessionInfo(), file = file.path("results", "session_info.txt"))
  
  # End timing
  end_time <- Sys.time()
  message("\nAnalysis completed at: ", end_time)
  message("Total runtime: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  
  return(seurat_obj)
}

# Configuration
config <- list(
  project_name = "FRC_scRNAseq_Analysis",
  species = "mouse",
  sample_dirs = c(
    "GFP-Ag_2D_LP" = "Data/aTomei_zWilkes_scRNAseq_12162021/multi_sample/GFP-Ag_2D_LP_sample_feature_bc_matrix",
    "GFP-Ag_3D_LP" = "Data/aTomei_zWilkes_scRNAseq_12162021/multi_sample/GFP-Ag_3D_LP_sample_feature_bc_matrix",
    "Fresh_FRCs" = "Data/aTomei_zWilkes_scRNAseq_12162021/single_sample/10X1028GeXSI_12162021_filtered_feature_bc_matrix"
  ),
  qc = list(
    min_genes = 200,
    max_genes = 6000,
    max_mt = 20,
    min_counts = 1000,
    max_counts = 50000
  )
)

# Run the analysis
if (sys.nframe() == 0) {
  tryCatch({
    seurat_obj <- run_analysis(config)
    message("Analysis completed successfully!")
  }, error = function(e) {
    message("Error in analysis: ", conditionMessage(e))
    traceback()
    quit(status = 1)
  }, finally = {
    # Close log file
    sink(type = "message")
    close(log_file)
  })
}
