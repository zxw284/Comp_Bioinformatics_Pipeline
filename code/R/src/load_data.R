#!/usr/bin/env Rscript
# Data Loading Module
# Description: Handles loading and initial processing of scRNA-seq data

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(AnnotationHub)
  library(ensembldb)
})

# Function to load 10X data and create Seurat objects
load_10x_data <- function(sample_dirs, sample_names, min_cells = 3, min_features = 200) {
  if (length(sample_dirs) != length(sample_names)) {
    stop("sample_dirs and sample_names must have the same length")
  }
  
  seurat_objects <- list()
  
  for (i in seq_along(sample_dirs)) {
    message("Loading sample: ", sample_names[i])
    
    # Check if the directory exists
    if (!dir.exists(sample_dirs[i])) {
      stop("Directory not found: ", sample_dirs[i])
    }
    
    # Load 10X data
    data <- Read10X(data.dir = sample_dirs[i])
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
      counts = if (is.list(data)) data$`Gene Expression` else data,
      project = sample_names[i],
      min.cells = min_cells,
      min.features = min_features
    )
    
    # Add sample metadata
    seurat_obj$sample_id <- sample_names[i]
    
    seurat_objects[[sample_names[i]]] <- seurat_obj
  }
  
  return(seurat_objects)
}

# Function to load cell cycle genes
get_cell_cycle_genes <- function(species = "Mus musculus") {
  # Try to load from local file first
  local_cc_file <- "data/cell_cycle_genes.csv"
  
  if (file.exists(local_cc_file)) {
    message("Loading cell cycle genes from local file")
    cc_genes <- read.csv(local_cc_file, stringsAsFactors = FALSE)
  } else {
    message("Downloading cell cycle genes from GitHub")
    cc_url <- "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv"
    cc_genes <- read.csv(cc_url, stringsAsFactors = FALSE)
    
    # Save locally for future use
    if (!dir.exists(dirname(local_cc_file))) {
      dir.create(dirname(local_cc_file), recursive = TRUE)
    }
    write.csv(cc_genes, local_cc_file, row.names = FALSE)
  }
  
  return(cc_genes)
}

# Function to get Ensembl gene annotations
get_gene_annotations <- function(species = "Mus musculus") {
  # Connect to AnnotationHub
  ah <- AnnotationHub()
  
  # Query for Ensembl database
  ahDb <- query(ah, pattern = c(species, "EnsDb"), ignore.case = TRUE)
  
  # Get the latest Ensembl database
  id <- ahDb %>% 
    mcols() %>% 
    rownames() %>% 
    tail(n = 1)
  
  # Load the database
  edb <- ah[[id]]
  
  # Extract gene-level information
  gene_annot <- genes(edb, return.type = "data.frame")
  
  return(gene_annot)
}

# Function to merge multiple Seurat objects
merge_seurat_objects <- function(seurat_objects, add_cell_ids = NULL, project = "merged") {
  if (is.null(add_cell_ids)) {
    add_cell_ids <- names(seurat_objects)
  }
  
  if (length(seurat_objects) == 1) {
    return(seurat_objects[[1]])
  }
  
  # Merge all objects
  merged_obj <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[-1],
    add.cell.ids = add_cell_ids,
    project = project
  )
  
  return(merged_obj)
}

# Main function to load all data
load_scrnaseq_data <- function(config) {
  # Create results directory if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
  }
  
  # Load cell cycle genes
  cell_cycle_genes <- get_cell_cycle_genes()
  
  # Get gene annotations
  gene_annotations <- get_gene_annotations()
  
  # Load 10X data
  seurat_objects <- load_10x_data(
    sample_dirs = config$sample_dirs,
    sample_names = names(config$sample_dirs)
  )
  
  # Merge if multiple samples
  if (length(seurat_objects) > 1) {
    seurat_obj <- merge_seurat_objects(
      seurat_objects,
      add_cell_ids = names(seurat_objects),
      project = config$project_name
    )
  } else {
    seurat_obj <- seurat_objects[[1]]
  }
  
  # Add metadata
  seurat_obj$orig.ident <- factor(seurat_obj$orig.ident)
  
  # Save the raw data
  saveRDS(seurat_obj, file = file.path("results", "seurat_raw.rds"))
  
  return(list(
    seurat = seurat_obj,
    cell_cycle_genes = cell_cycle_genes,
    gene_annotations = gene_annotations
  ))
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Example configuration
  config <- list(
    sample_dirs = c(
      "GFP-Ag_2D_LP" = "Data/aTomei_zWilkes_scRNAseq_12162021/multi_sample/GFP-Ag_2D_LP_sample_feature_bc_matrix",
      "GFP-Ag_3D_LP" = "Data/aTomei_zWilkes_scRNAseq_12162021/multi_sample/GFP-Ag_3D_LP_sample_feature_bc_matrix",
      "Fresh_FRCs" = "Data/aTomei_zWilkes_scRNAseq_12162021/single_sample/10X1028GeXSI_12162021_filtered_feature_bc_matrix"
    ),
    project_name = "FRC_scRNAseq_Analysis"
  )
  
  # Load data
  data <- load_scrnaseq_data(config)
  print("Data loading complete!")
}
