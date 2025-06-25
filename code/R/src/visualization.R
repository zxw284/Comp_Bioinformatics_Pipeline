#!/usr/bin/env Rscript
# Visualization Module
# Description: Provides visualization functions for scRNA-seq analysis

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
  library(ggridges)
  library(viridis)
  library(pheatmap)
  library(ggpubr)
})

# Set default ggplot2 theme
set_default_theme <- function() {
  theme_set(
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )
  )
}

# Create a color palette
get_palette <- function(palette = "Set1", n = 10, alpha = 1) {
  if (palette %in% rownames(brewer.pal.info)) {
    max_colors <- brewer.pal.info[palette, "maxcolors"]
    if (n <= max_colors) {
      cols <- brewer.pal(n, palette)
    } else {
      cols <- colorRampPalette(brewer.pal(max_colors, palette))(n)
    }
  } else {
    cols <- viridis(n, option = palette, alpha = alpha)
  }
  return(cols)
}

# Plot UMAP/tSNE with custom aesthetics
plot_dim_reduction <- function(seurat_obj,
                             reduction = "umap",
                             group.by = "seurat_clusters",
                             label = TRUE,
                             label.size = 4,
                             pt.size = 0.5,
                             alpha = 0.7,
                             palette = "Set1",
                             legend.position = "right",
                             legend.title = NULL,
                             title = NULL) {
  
  # Get reduction coordinates
  if (is.null(seurat_obj@reductions[[reduction]])) {
    stop(paste("Reduction", reduction, "not found in the Seurat object"))
  }
  
  # Create data frame for plotting
  plot_data <- FetchData(seurat_obj, vars = c(
    paste0(reduction, "_1"),
    paste0(reduction, "_2"),
    group.by
  ))
  
  colnames(plot_data) <- c("x", "y", "group")
  
  # Get cluster centers for labels
  if (label) {
    centers <- plot_data %>%
      group_by(group) %>%
      summarize(
        x = median(x),
        y = median(y)
      )
  }
  
  # Create color palette
  n_groups <- length(unique(plot_data$group))
  colors <- get_palette(palette, n = n_groups)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_point(size = pt.size, alpha = alpha) +
    scale_color_manual(values = colors) +
    labs(
      x = paste(toupper(reduction), "1"),
      y = paste(toupper(reduction), "2"),
      color = if (is.null(legend.title)) group.by else legend.title,
      title = title
    ) +
    theme_void() +
    theme(
      legend.position = legend.position,
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Add cluster labels
  if (label) {
    p <- p + geom_text_repel(
      data = centers,
      aes(x = x, y = y, label = group),
      size = label.size,
      color = "black",
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.size = 0.2
    )
  }
  
  return(p)
}

# Plot feature expression on UMAP/tSNE
plot_feature_expression <- function(seurat_obj,
                                  features,
                                  reduction = "umap",
                                  slot = "data",
                                  ncol = NULL,
                                  pt.size = 0.5,
                                  alpha = 0.7,
                                  order = TRUE,
                                  min.cutoff = NA,
                                  max.cutoff = NA,
                                  title = NULL) {
  
  # Check if features exist
  features <- features[features %in% rownames(seurat_obj)]
  if (length(features) == 0) {
    stop("None of the specified features were found in the Seurat object")
  }
  
  # Get reduction coordinates
  if (is.null(seurat_obj@reductions[[reduction]])) {
    stop(paste("Reduction", reduction, "not found in the Seurat object"))
  }
  
  # Create plots
  plot_list <- lapply(features, function(feature) {
    # Get expression data
    expr <- FetchData(seurat_obj, vars = feature, slot = slot)[, 1]
    
    # Apply cutoffs
    if (!is.na(min.cutoff)) {
      expr[expr < min.cutoff] <- min.cutoff
    }
    if (!is.na(max.cutoff)) {
      expr[expr > max.cutoff] <- max.cutoff
    }
    
    # Create data frame for plotting
    plot_data <- data.frame(
      x = seurat_obj@reductions[[reduction]]@cell.embeddings[, 1],
      y = seurat_obj@reductions[[reduction]]@cell.embeddings[, 2],
      expression = expr
    )
    
    # Order points by expression if requested
    if (order) {
      plot_data <- plot_data[order(plot_data$expression, decreasing = FALSE), ]
    }
    
    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = y, color = expression)) +
      geom_point(size = pt.size, alpha = alpha) +
      scale_color_viridis(option = "C") +
      labs(
        x = paste(toupper(reduction), "1"),
        y = paste(toupper(reduction), "2"),
        title = if (!is.null(title)) title else feature,
        color = "Expression"
      ) +
      theme_void() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    return(p)
  })
  
  # Combine plots if multiple features
  if (length(plot_list) > 1) {
    if (is.null(ncol)) {
      ncol <- min(3, length(plot_list))
    }
    return(wrap_plots(plot_list, ncol = ncol))
  } else {
    return(plot_list[[1]])
  }
}

# Create a violin plot of gene expression
plot_violin <- function(seurat_obj,
                       features,
                       group.by = "seurat_clusters",
                       slot = "data",
                       ncol = NULL,
                       pt.size = 0,
                       alpha = 0.7,
                       add.boxplot = TRUE,
                       colors = NULL,
                       ylab = "Expression level") {
  
  # Check if features exist
  features <- features[features %in% rownames(seurat_obj)]
  if (length(features) == 0) {
    stop("None of the specified features were found in the Seurat object")
  }
  
  # Get expression data
  expr_data <- FetchData(seurat_obj, vars = c(features, group.by))
  
  # Melt data for plotting
  plot_data <- reshape2::melt(
    expr_data,
    id.vars = group.by,
    variable.name = "feature",
    value.name = "expression"
  )
  
  # Set colors if not provided
  n_groups <- length(unique(plot_data[[group.by]]))
  if (is.null(colors)) {
    colors <- get_palette("Set1", n = n_groups)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes_string(x = group.by, y = "expression", fill = group.by)) +
    geom_violin(scale = "width", alpha = alpha) +
    scale_fill_manual(values = colors) +
    facet_wrap(~ feature, scales = "free_y", ncol = ncol) +
    labs(x = "", y = ylab) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  # Add boxplot if requested
  if (add.boxplot) {
    p <- p + geom_boxplot(width = 0.1, fill = "white", alpha = 0.5)
  }
  
  # Add points if requested
  if (pt.size > 0) {
    p <- p + geom_jitter(width = 0.2, size = pt.size, alpha = 0.5)
  }
  
  return(p)
}

# Create a heatmap of marker genes
plot_marker_heatmap <- function(seurat_obj,
                              markers,
                              group.by = "seurat_clusters",
                              top_n = 10,
                              slot = "scale.data",
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              show_rownames = TRUE,
                              show_colnames = FALSE,
                              fontsize_row = 8,
                              fontsize_col = 8,
                              ...) {
  
  # Get top markers
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC)
  
  # Get expression matrix
  expr_mat <- GetAssayData(seurat_obj, slot = slot)
  
  # Subset to marker genes
  expr_mat <- expr_mat[unique(top_markers$gene), ]
  
  # Get cluster annotations
  annotation_col <- data.frame(
    Cluster = seurat_obj@meta.data[[group.by]]
  )
  rownames(annotation_col) <- colnames(seurat_obj)
  
  # Create color palette
  cluster_colors <- get_palette("Set1", n = length(unique(annotation_col$Cluster)))
  names(cluster_colors) <- unique(annotation_col$Cluster)
  
  # Create annotation colors
  ann_colors <- list(
    Cluster = cluster_colors
  )
  
  # Create heatmap
  pheatmap(
    mat = expr_mat,
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
    scale = "row",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    ...
  )
}

# Create a dot plot of marker genes
plot_dot_plot <- function(seurat_obj,
                        markers,
                        group.by = "seurat_clusters",
                        top_n = 3,
                        cols = c("lightgrey", "blue"),
                        dot.scale = 6,
                        ...) {
  
  # Get top markers
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC)
  
  # Create dot plot
  DotPlot(
    object = seurat_obj,
    features = unique(top_markers$gene),
    group.by = group.by,
    cols = cols,
    dot.scale = dot.scale,
    ...
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_blank()
    )
}

# Save plot with consistent styling
save_plot <- function(plot,
                     filename,
                     width = 10,
                     height = 8,
                     dpi = 300,
                     ...) {
  
  # Create directory if it doesn't exist
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  # Get file extension
  ext <- tools::file_ext(filename)
  
  # Save plot
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  
  message("Plot saved to: ", filename)
}

# Set default theme when the package is loaded
set_default_theme()

# Run if executed directly
if (sys.nframe() == 0) {
  # Example usage
  # library(Seurat)
  # 
  # # Load data
  # seurat_obj <- readRDS("results/clustering/seurat_clustered.rds")
  # 
  # # Create UMAP plot
  # p1 <- plot_dim_reduction(
  #   seurat_obj = seurat_obj,
  #   reduction = "umap",
  #   group.by = "seurat_clusters",
  #   title = "UMAP by Cluster"
  # )
  # 
  # # Save plot
  # save_plot(
  #   plot = p1,
  #   filename = "results/plots/umap_clusters.png",
  #   width = 10,
  #   height = 8
  # )
}
