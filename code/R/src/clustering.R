#!/usr/bin/env Rscript
# Clustering and Cell Type Annotation Module
# Description: Performs clustering and cell type annotation on scRNA-seq data

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(clustree)
  library(SingleR)
  library(celldex)
  library(scater)
})

# Function to find neighbors
find_neighbors <- function(seurat_obj, 
                         reduction = "pca",
                         dims = 1:30,
                         k.param = 20,
                         prune.SNN = 1/15,
                         force.recalc = FALSE) {
  
  seurat_obj <- FindNeighbors(
    object = seurat_obj,
    reduction = reduction,
    dims = dims,
    k.param = k.param,
    prune.SNN = prune.SNN,
    force.recalc = force.recalc,
    verbose = FALSE
  )
  
  return(seurat_obj)
}

# Function to find clusters
find_clusters <- function(seurat_obj,
                        resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),
                        algorithm = 1,
                        random.seed = 42) {
  
  # Set seed for reproducibility
  set.seed(random.seed)
  
  # Find clusters at multiple resolutions
  seurat_obj <- FindClusters(
    object = seurat_obj,
    resolution = resolution,
    algorithm = algorithm,
    verbose = FALSE
  )
  
  return(seurat_obj)
}

# Function to visualize cluster stability
visualize_cluster_stability <- function(seurat_obj, 
                                      prefix = "SCT_snn_res.",
                                      output_dir = "results/clustering") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get clustering results
  clustering_columns <- grep(paste0("^", prefix), colnames(seurat_obj@meta.data), value = TRUE)
  
  if (length(clustering_columns) < 2) {
    warning("Need at least two clustering resolutions to compare.")
    return(NULL)
  }
  
  # Create clustree plot
  clustree_plot <- clustree(seurat_obj, prefix = prefix) +
    theme(legend.position = "bottom") +
    ggtitle("Cluster Stability Across Resolutions")
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, "clustree_plot.png"),
    plot = clustree_plot,
    width = 10,
    height = 12,
    dpi = 300
  )
  
  return(clustree_plot)
}

# Function to annotate cell types using SingleR
annotate_cell_types <- function(seurat_obj, 
                              ref = c("HumanPrimaryCellAtlasData", "DatabaseImmuneCellExpressionData", "NovershternHematopoieticData"),
                              clusters = "seurat_clusters",
                              labels = NULL,
                              species = "mouse") {
  
  # Convert species to match celldex format
  if (tolower(species) == "mouse") {
    ref <- sub("Data$", "Data('mouse')", ref)
  } else if (tolower(species) == "human") {
    ref <- sub("Data$", "Data('human')", ref)
  }
  
  # Load reference data
  ref_data <- lapply(ref, function(r) eval(parse(text = r)))
  names(ref_data) <- ref
  
  # Get expression data
  expr <- GetAssayData(seurat_obj, slot = "data")
  
  # Get cluster labels if not provided
  if (is.null(labels)) {
    if (clusters %in% colnames(seurat_obj@meta.data)) {
      labels <- seurat_obj@meta.data[[clusters]]
    } else {
      stop("Cluster column not found in metadata. Please specify a valid column name.")
    }
  }
  
  # Run SingleR for each reference
  pred_list <- list()
  for (i in seq_along(ref_data)) {
    message("Running SingleR with reference: ", names(ref_data)[i])
    
    # Convert to common genes
    common_genes <- intersect(rownames(expr), rownames(ref_data[[i]]))
    
    # Run SingleR
    pred <- SingleR(
      test = expr[common_genes, ],
      ref = ref_data[[i]][common_genes, ],
      labels = ref_data[[i]]$label.fine,
      clusters = labels,
      de.method = "wilcox"
    )
    
    pred_list[[names(ref_data)[i]]] <- pred
  }
  
  # Get consensus annotations
  consensus_labels <- lapply(pred_list, function(x) x$labels)
  consensus_df <- do.call(cbind, consensus_labels)
  
  # Get most frequent label for each cell
  final_labels <- apply(consensus_df, 1, function(x) {
    tab <- table(x)
    names(tab)[which.max(tab)]
  })
  
  # Add to metadata
  seurat_obj$cell_type <- final_labels[as.character(seurat_obj@meta.data[[clusters]])]
  
  # Return results
  return(list(
    seurat = seurat_obj,
    predictions = pred_list
  ))
}

# Function to visualize cell type proportions
visualize_cell_type_proportions <- function(seurat_obj,
                                          group_by = "orig.ident",
                                          cell_type_col = "cell_type",
                                          output_dir = "results/clustering") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Calculate proportions
  prop_df <- seurat_obj@meta.data %>%
    group_by(!!sym(group_by), !!sym(cell_type_col)) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(!!sym(group_by)) %>%
    mutate(proportion = count / sum(count) * 100)
  
  # Create stacked bar plot
  prop_plot <- ggplot(prop_df, aes(x = !!sym(group_by), y = proportion, fill = !!sym(cell_type_col))) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "Sample",
      y = "Proportion (%)",
      fill = "Cell Type",
      title = "Cell Type Proportions by Sample"
    )
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, "cell_type_proportions.png"),
    plot = prop_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  return(prop_plot)
}

# Main function for clustering and annotation
perform_clustering <- function(seurat_obj,
                             reduction = "pca",
                             dims = 1:30,
                             resolution = seq(0.2, 2.0, by = 0.2),
                             k.param = 20,
                             algorithm = 1,
                             annotate = TRUE,
                             ref_datasets = c("HumanPrimaryCellAtlasData", "DatabaseImmuneCellExpressionData"),
                             species = "mouse",
                             output_dir = "results/clustering") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Finding nearest neighbors...")
  seurat_obj <- find_neighbors(
    seurat_obj = seurat_obj,
    reduction = reduction,
    dims = dims,
    k.param = k.param
  )
  
  message("Finding clusters...")
  seurat_obj <- find_clusters(
    seurat_obj = seurat_obj,
    resolution = resolution,
    algorithm = algorithm
  )
  
  # Visualize cluster stability
  message("Visualizing cluster stability...")
  clustree_plot <- visualize_cluster_stability(
    seurat_obj = seurat_obj,
    output_dir = output_dir
  )
  
  # Set default resolution if not specified
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    default_res <- 0.8  # Default resolution
    message(sprintf("No default clustering found. Using resolution = %.1f", default_res))
    seurat_obj$seurat_clusters <- seurat_obj@meta.data[[sprintf("SCT_snn_res.%.1f", default_res)]]
  }
  
  # Annotate cell types if requested
  if (annotate) {
    message("Annotating cell types...")
    annotation_result <- annotate_cell_types(
      seurat_obj = seurat_obj,
      ref = ref_datasets,
      species = species
    )
    
    seurat_obj <- annotation_result$seurat
    
    # Visualize cell type proportions
    message("Visualizing cell type proportions...")
    prop_plot <- visualize_cell_type_proportions(
      seurat_obj = seurat_obj,
      output_dir = output_dir
    )
  }
  
  # Save results
  saveRDS(seurat_obj, file.path(output_dir, "seurat_clustered.rds"))
  
  # Return results
  result <- list(
    seurat = seurat_obj,
    clustree_plot = if (exists("clustree_plot")) clustree_plot else NULL,
    prop_plot = if (exists("prop_plot")) prop_plot else NULL
  )
  
  if (exists("annotation_result")) {
    result$annotation_result <- annotation_result$predictions
  }
  
  return(result)
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Example usage
  # library(Seurat)
  # 
  # # Load data
  # seurat_obj <- readRDS("results/dimensionality_reduction/seurat_dim_reduced.rds")
  # 
  # # Run clustering and annotation
  # clustering_result <- perform_clustering(
  #   seurat_obj = seurat_obj,
  #   dims = 1:30,
  #   resolution = seq(0.2, 2.0, by = 0.2),
  #   annotate = TRUE,
  #   species = "mouse"
  # )
  # 
  # print("Clustering and annotation complete!")
}
