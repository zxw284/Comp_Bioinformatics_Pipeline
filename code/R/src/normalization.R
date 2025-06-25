#!/usr/bin/env Rscript
# Normalization Module
# Description: Handles data normalization and feature selection for scRNA-seq data

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)
  library(future)
  library(dplyr)
})

# Function to normalize data using SCTransform
normalize_sctransform <- function(seurat_obj, 
                                 vars_to_regress = NULL, 
                                 verbose = TRUE,
                                 ncells = 5000,
                                 n_genes = 3000,
                                 conserve.memory = FALSE) {
  
  # Set up parallel processing if available
  if (future::supportsMulticore()) {
    future::plan("multicore", workers = future::availableCores() - 1)
  } else {
    future::plan("multisession")
  }
  
  # Set verbosity
  options(future.globals.maxSize = 8000 * 1024^2)  # 8GB per worker
  
  # Run SCTransform
  seurat_obj <- SCTransform(
    object = seurat_obj,
    method = "glmGamPoi",
    vars.to.regress = vars_to_regress,
    verbose = verbose,
    ncells = ncells,
    n_genes = n_genes,
    conserve.memory = conserve.memory,
    return.only.var.genes = TRUE,
    variable.features.n = 3000,
    vst.flavor = "v2"
  )
  
  return(seurat_obj)
}

# Function to select variable features
select_variable_features <- function(seurat_obj, 
                                   nfeatures = 2000, 
                                   selection.method = "vst") {
  
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj,
    selection.method = selection.method,
    nfeatures = nfeatures,
    verbose = FALSE
  )
  
  return(seurat_obj)
}

# Function to scale data
scale_data <- function(seurat_obj, 
                      features = NULL, 
                      vars.to.regress = NULL,
                      verbose = TRUE) {
  
  seurat_obj <- ScaleData(
    object = seurat_obj,
    features = features,
    vars.to.regress = vars.to.regress,
    verbose = verbose
  )
  
  return(seurat_obj)
}

# Function to regress out cell cycle effects
regress_cell_cycle <- function(seurat_obj, 
                              s.features, 
                              g2m.features,
                              species = "mouse") {
  
  # Convert gene names to appropriate case
  if (species == "mouse") {
    s.features <- stringr::str_to_title(s.features)
    g2m.features <- stringr::str_to_title(g2m.features)
  }
  
  # Score cell cycle
  seurat_obj <- CellCycleScoring(
    object = seurat_obj,
    s.features = s.features,
    g2m.features = g2m.features,
    set.ident = FALSE
  )
  
  # Regress out cell cycle scores
  seurat_obj <- ScaleData(
    object = seurat_obj,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(seurat_obj)
  )
  
  return(seurat_obj)
}

# Main function for data normalization
normalize_data <- function(seurat_obj, 
                         method = "sctransform",
                         vars_to_regress = NULL,
                         cell_cycle_genes = NULL,
                         species = "mouse",
                         output_dir = "results/normalization") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Apply normalization
  if (method == "sctransform") {
    message("Performing SCTransform normalization...")
    seurat_obj <- normalize_sctransform(
      seurat_obj = seurat_obj,
      vars_to_regress = vars_to_regress
    )
  } else {
    stop("Only 'sctransform' method is currently supported")
  }
  
  # Regress out cell cycle effects if genes are provided
  if (!is.null(cell_cycle_genes)) {
    message("Regressing out cell cycle effects...")
    seurat_obj <- regress_cell_cycle(
      seurat_obj = seurat_obj,
      s.features = cell_cycle_genes$s.genes,
      g2m.features = cell_cycle_genes$g2m.genes,
      species = species
    )
  }
  
  # Save variable features plot
  var_feat_plot <- VariableFeaturePlot(seurat_obj) + 
    theme(legend.position = "bottom") +
    ggtitle("Variable Features")
  
  ggsave(
    filename = file.path(output_dir, "variable_features.png"),
    plot = var_feat_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Save normalized object
  saveRDS(seurat_obj, file.path(output_dir, "seurat_normalized.rds"))
  
  return(seurat_obj)
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Example usage
  # library(Seurat)
  # 
  # # Load data
  # seurat_obj <- readRDS("results/QC/seurat_qc_filtered.rds")
  # 
  # # Load cell cycle genes if needed
  # # cell_cycle_genes <- readRDS("data/cell_cycle_genes.rds")
  # 
  # # Normalize data
  # seurat_norm <- normalize_data(
  #   seurat_obj = seurat_obj,
  #   method = "sctransform",
  #   vars_to_regress = c("percent.mt"),
  #   cell_cycle_genes = NULL,
  #   species = "mouse"
  # )
  # 
  # print("Data normalization complete!")
}
