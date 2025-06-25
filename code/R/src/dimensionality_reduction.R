#!/usr/bin/env Rscript
# Dimensionality Reduction Module
# Description: Performs PCA, t-SNE, and UMAP on scRNA-seq data

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(umap)
})

# Function to run PCA
run_pca <- function(seurat_obj, 
                   npcs = 50, 
                   features = NULL,
                   verbose = TRUE) {
  
  # If features not specified, use variable features
  if (is.null(features)) {
    features <- VariableFeatures(seurat_obj)
    if (length(features) == 0) {
      stop("No variable features found. Run FindVariableFeatures() first.")
    }
  }
  
  # Run PCA
  seurat_obj <- RunPCA(
    object = seurat_obj,
    features = features,
    npcs = npcs,
    verbose = verbose
  )
  
  return(seurat_obj)
}

# Function to determine significant PCs
choose_pcs <- function(seurat_obj, 
                      reduction = "pca",
                      method = "elbow",
                      ndims = 50) {
  
  # Elbow plot method
  if (method == "elbow") {
    pct_var <- seurat_obj[[reduction]]@stdev^2 / sum(seurat_obj[[reduction]]@stdev^2) * 100
    cum_var <- cumsum(pct_var)
    
    # Find the elbow point
    d1 <- diff(pct_var)
    d2 <- diff(d1)
    elbow <- which.min(d2[1:(length(d2)-1)] / d1[2:(length(d2))]) + 1
    
    # Create elbow plot
    elbow_plot <- data.frame(
      PC = 1:length(pct_var),
      Variance = pct_var,
      Cumulative = cum_var
    )
    
    p1 <- ggplot(elbow_plot, aes(x = PC, y = Variance)) +
      geom_point() +
      geom_vline(xintercept = elbow, linetype = "dashed", color = "red") +
      ggtitle("Elbow Plot") +
      theme_bw()
    
    p2 <- ggplot(elbow_plot, aes(x = PC, y = Cumulative)) +
      geom_point() +
      geom_hline(yintercept = 90, linetype = "dashed", color = "blue") +
      ylim(0, 100) +
      ggtitle("Cumulative Variance") +
      theme_bw()
    
    # Combine plots
    pc_plot <- p1 / p2 + plot_annotation(
      title = "Principal Component Analysis",
      subtitle = paste0("Suggested number of PCs: ", elbow)
    )
    
    return(list(
      n_pcs = elbow,
      plot = pc_plot,
      pct_var = pct_var,
      cum_var = cum_var
    ))
    
  } else if (method == "jackstraw") {
    # JackStraw method
    seurat_obj <- JackStraw(seurat_obj, num.replicate = 100, dims = ndims)
    seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:ndims)
    
    # Return JackStraw plot and significant PCs
    return(list(
      plot = JackStrawPlot(seurat_obj, dims = 1:ndims),
      scores = seurat_obj@reductions$pca@jackstraw$overall.p.values
    ))
  }
}

# Function to run t-SNE
run_tsne <- function(seurat_obj, 
                    dims = 1:30, 
                    reduction = "pca",
                    perplexity = 30,
                    seed = 42) {
  
  set.seed(seed)
  seurat_obj <- RunTSNE(
    object = seurat_obj,
    reduction = reduction,
    dims = dims,
    perplexity = perplexity,
    verbose = FALSE
  )
  
  return(seurat_obj)
}

# Function to run UMAP
run_umap <- function(seurat_obj, 
                    dims = 1:30, 
                    reduction = "pca",
                    n.neighbors = 30,
                    min.dist = 0.3,
                    metric = "cosine",
                    seed = 42) {
  
  set.seed(seed)
  seurat_obj <- RunUMAP(
    object = seurat_obj,
    reduction = reduction,
    dims = dims,
    n.neighbors = n.neighbors,
    min.dist = min.dist,
    metric = metric,
    verbose = FALSE
  )
  
  return(seurat_obj)
}

# Function to visualize embeddings
plot_embeddings <- function(seurat_obj, 
                           reduction = "umap",
                           group_by = NULL,
                           split_by = NULL,
                           label = FALSE,
                           colors = NULL,
                           pt.size = 0.5,
                           alpha = 0.7) {
  
  # Get available reductions
  available_reductions <- names(seurat_obj@reductions)
  if (!reduction %in% available_reductions) {
    stop("Specified reduction not found. Available reductions: ", 
         paste(available_reductions, collapse = ", "))
  }
  
  # Create base plot
  p <- DimPlot(
    object = seurat_obj,
    reduction = reduction,
    group.by = group_by,
    split.by = split_by,
    label = label,
    pt.size = pt.size,
    label.size = 3,
    repel = TRUE
  ) + 
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  # Adjust transparency
  if (alpha != 1) {
    p$layers[[1]]$aes_params$alpha <- alpha
  }
  
  return(p)
}

# Main function for dimensionality reduction
perform_dimensionality_reduction <- function(seurat_obj,
                                           npcs = 50,
                                           pca_features = NULL,
                                           reduction_methods = c("pca", "umap"),
                                           tsne_perplexity = 30,
                                           umap_n_neighbors = 30,
                                           umap_min_dist = 0.3,
                                           group_by = "orig.ident",
                                           output_dir = "results/dimensionality_reduction") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run PCA if requested
  if (any(c("pca", "tsne", "umap") %in% tolower(reduction_methods))) {
    message("Running PCA...")
    seurat_obj <- run_pca(
      seurat_obj = seurat_obj,
      npcs = npcs,
      features = pca_features
    )
    
    # Plot PCA results
    pca_plots <- list()
    pca_plots$elbow <- ElbowPlot(seurat_obj, ndims = npcs) + 
      ggtitle("Elbow Plot") + 
      theme_bw()
    
    pca_plots$pca_heatmap <- DimHeatmap(
      seurat_obj, 
      dims = 1:min(15, npcs), 
      cells = 500, 
      balanced = TRUE,
      fast = FALSE
    )
    
    # Save PCA plots
    ggsave(
      filename = file.path(output_dir, "pca_elbow_plot.png"),
      plot = pca_plots$elbow,
      width = 8,
      height = 6,
      dpi = 300
    )
    
    # Determine significant PCs
    pc_stats <- choose_pcs(seurat_obj, method = "elbow", ndims = npcs)
    n_pcs <- pc_stats$n_pcs
    
    message("Using ", n_pcs, " principal components for downstream analysis.")
  }
  
  # Run t-SNE if requested
  if ("tsne" %in% tolower(reduction_methods)) {
    message("Running t-SNE...")
    seurat_obj <- run_tsne(
      seurat_obj = seurat_obj,
      dims = 1:n_pcs,
      perplexity = tsne_perplexity
    )
    
    # Plot t-SNE
    tsne_plot <- plot_embeddings(
      seurat_obj = seurat_obj,
      reduction = "tsne",
      group_by = group_by,
      label = TRUE
    )
    
    # Save t-SNE plot
    ggsave(
      filename = file.path(output_dir, "tsne_plot.png"),
      plot = tsne_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  # Run UMAP if requested
  if ("umap" %in% tolower(reduction_methods)) {
    message("Running UMAP...")
    seurat_obj <- run_umap(
      seurat_obj = seurat_obj,
      dims = 1:n_pcs,
      n.neighbors = umap_n_neighbors,
      min.dist = umap_min_dist
    )
    
    # Plot UMAP
    umap_plot <- plot_embeddings(
      seurat_obj = seurat_obj,
      reduction = "umap",
      group_by = group_by,
      label = TRUE
    )
    
    # Save UMAP plot
    ggsave(
      filename = file.path(output_dir, "umap_plot.png"),
      plot = umap_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  # Save results
  saveRDS(seurat_obj, file.path(output_dir, "seurat_dim_reduced.rds"))
  
  # Return results
  return(seurat_obj)
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Example usage
  # library(Seurat)
  # 
  # # Load data
  # seurat_obj <- readRDS("results/normalization/seurat_normalized.rds")
  # 
  # # Run dimensionality reduction
  # seurat_dim_red <- perform_dimensionality_reduction(
  #   seurat_obj = seurat_obj,
  #   npcs = 50,
  #   reduction_methods = c("pca", "tsne", "umap"),
  #   group_by = "orig.ident"
  # )
  # 
  # print("Dimensionality reduction complete!")
}
