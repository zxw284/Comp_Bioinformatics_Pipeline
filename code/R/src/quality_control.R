#!/usr/bin/env Rscript
# Quality Control Module
# Description: Performs quality control on single-cell RNA-seq data

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
})

# Calculate mitochondrial and ribosomal percentages
calculate_qc_metrics <- function(seurat_obj, species = "mouse") {
  # Set species-specific patterns
  if (tolower(species) == "mouse") {
    mt_pattern <- "^mt-"
    ribo_pattern <- "^R[Pp][LlSs]"
  } else {
    mt_pattern <- "^MT-"
    ribo_pattern <- "^R[Pp][LlSs]"
  }
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  # Calculate ribosomal percentage
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = ribo_pattern)
  
  # Calculate log10 genes per UMI
  seurat_obj[["log10_genes_per_UMI"]] <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  return(seurat_obj)
}

# Generate QC plots
plot_qc_metrics <- function(seurat_obj, output_dir = "results/QC") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Violin plots of QC metrics
  vln_plot <- VlnPlot(seurat_obj, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                     ncol = 2, 
                     pt.size = 0.1) &
    theme(plot.title = element_text(size = 10))
  
  # Scatter plots
  scatter1 <- FeatureScatter(seurat_obj, 
                           feature1 = "nCount_RNA", 
                           feature2 = "nFeature_RNA",
                           pt.size = 0.5) +
    theme(legend.position = "none")
  
  scatter2 <- FeatureScatter(seurat_obj, 
                           feature1 = "nCount_RNA", 
                           feature2 = "percent.mt",
                           pt.size = 0.5) +
    theme(legend.position = "none")
  
  scatter3 <- FeatureScatter(seurat_obj, 
                           feature1 = "nFeature_RNA", 
                           feature2 = "percent.mt",
                           pt.size = 0.5)
  
  # Combine plots
  combined_plot <- (vln_plot / (scatter1 + scatter2) / scatter3) + 
    plot_annotation(title = "QC Metrics")
  
  # Save plots
  ggsave(file.path(output_dir, "qc_metrics.png"), 
         combined_plot, 
         width = 12, 
         height = 16,
         dpi = 300)
  
  # Save individual plots
  ggsave(file.path(output_dir, "qc_violin_plots.png"), 
         vln_plot, 
         width = 10, 
         height = 8,
         dpi = 300)
  
  ggsave(file.path(output_dir, "qc_scatter_plots.png"), 
         scatter1 + scatter2 + scatter3, 
         width = 15, 
         height = 5,
         dpi = 300)
  
  return(combined_plot)
}

# Filter cells based on QC metrics
filter_cells <- function(seurat_obj, 
                        min_genes = 200, 
                        max_genes = 6000,
                        max_mt = 20,
                        min_counts = 1000,
                        max_counts = 50000) {
  
  message("Number of cells before filtering: ", ncol(seurat_obj))
  
  # Store filtering parameters in metadata
  seurat_obj@misc$filtering_params <- list(
    min_genes = min_genes,
    max_genes = max_genes,
    max_mt = max_mt,
    min_counts = min_counts,
    max_counts = max_counts,
    date = Sys.time()
  )
  
  # Apply filters
  cells_keep <- (seurat_obj$nFeature_RNA >= min_genes) &
    (seurat_obj$nFeature_RNA <= max_genes) &
    (seurat_obj$percent.mt <= max_mt) &
    (seurat_obj$nCount_RNA >= min_counts) &
    (seurat_obj$nCount_RNA <= max_counts)
  
  seurat_obj <- seurat_obj[, cells_keep]
  
  message("Number of cells after filtering: ", ncol(seurat_obj))
  message("Percentage of cells kept: ", 
          round(ncol(seurat_obj) / sum(cells_keep) * 100, 2), "%")
  
  return(seurat_obj)
}

# Run cell cycle scoring
score_cell_cycle <- function(seurat_obj, s_genes, g2m_genes, species = "mouse") {
  # Convert gene names to appropriate case
  if (species == "mouse") {
    s_genes <- str_to_title(s_genes)
    g2m_genes <- str_to_title(g2m_genes)
  }
  
  # Score cell cycle phases
  seurat_obj <- CellCycleScoring(
    object = seurat_obj,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = FALSE
  )
  
  return(seurat_obj)
}

# Main QC function
perform_quality_control <- function(seurat_obj, 
                                  cell_cycle_genes = NULL,
                                  output_dir = "results/QC",
                                  ...) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Calculating QC metrics...")
  seurat_obj <- calculate_qc_metrics(seurat_obj, ...)
  
  message("Generating QC plots...")
  qc_plots <- plot_qc_metrics(seurat_obj, output_dir)
  
  message("Filtering cells...")
  seurat_obj <- filter_cells(seurat_obj, ...)
  
  # Score cell cycle if genes are provided
  if (!is.null(cell_cycle_genes)) {
    message("Scoring cell cycle...")
    seurat_obj <- score_cell_cycle(
      seurat_obj,
      s_genes = cell_cycle_genes$s.genes,
      g2m_genes = cell_cycle_genes$g2m.genes,
      ...
    )
  }
  
  # Save filtered object
  saveRDS(seurat_obj, file.path(output_dir, "seurat_qc_filtered.rds"))
  
  return(seurat_obj)
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Example usage
  library(Seurat)
  
  # Load data (example)
  # seurat_obj <- readRDS("results/seurat_raw.rds")
  
  # Load cell cycle genes
  # cell_cycle_genes <- readRDS("data/cell_cycle_genes.rds")
  
  # Run QC
  # seurat_qc <- perform_quality_control(
  #   seurat_obj,
  #   cell_cycle_genes = cell_cycle_genes,
  #   species = "mouse",
  #   min_genes = 200,
  #   max_genes = 6000,
  #   max_mt = 20,
  #   min_counts = 1000,
  #   max_counts = 50000
  # )
  
  message("Quality control complete!")
}
