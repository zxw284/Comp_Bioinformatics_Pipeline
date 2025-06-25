#!/usr/bin/env Rscript
# Differential Expression Analysis Module
# Description: Performs differential expression analysis on scRNA-seq data

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(MAST)
  library(limma)
  library(DESeq2)
  library(EnhancedVolcano)
})

# Function to find differentially expressed genes between groups
find_markers <- function(seurat_obj,
                       ident.1 = NULL,
                       ident.2 = NULL,
                       group.by = NULL,
                       test.use = "wilcox",
                       min.pct = 0.1,
                       logfc.threshold = 0.25,
                       only.pos = FALSE,
                       max.cells.per.ident = Inf,
                       random.seed = 42) {
  
  # Set seed for reproducibility
  set.seed(random.seed)
  
  # Set identity if group.by is specified
  if (!is.null(group.by)) {
    Idents(seurat_obj) <- group.by
  }
  
  # Find markers
  markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos,
    max.cells.per.ident = max.cells.per.ident,
    verbose = FALSE
  )
  
  # Add gene names as a column
  markers$gene <- rownames(markers)
  
  # Calculate average expression
  if (!is.null(ident.1) && !is.null(ident.2)) {
    avg_expr <- AverageExpression(
      seurat_obj,
      features = markers$gene,
      group.by = group.by,
      slot = "data"
    )[[1]]
    
    # Add average expression to results
    markers$avg_expr1 <- rowMeans(avg_expr[markers$gene, ident.1, drop = FALSE])
    markers$avg_expr2 <- rowMeans(avg_expr[markers$gene, ident.2, drop = FALSE])
  }
  
  # Add significance
  markers$significant <- ifelse(markers$p_val_adj < 0.05, "Yes", "No")
  
  return(markers)
}

# Function to create a volcano plot
create_volcano_plot <- function(markers,
                              title = "Volcano Plot",
                              pval_cutoff = 0.05,
                              lfc_cutoff = 0.25,
                              top_n = 10,
                              label_all = FALSE) {
  
  # Prepare data
  markers$log_pval <- -log10(markers$p_val_adj)
  markers$log_pval[is.infinite(markers$log_pval)] <- max(markers$log_pval[is.finite(markers$log_pval)]) + 5
  
  # Highlight top genes
  markers$highlight <- "Not significant"
  markers$highlight[markers$p_val_adj < pval_cutoff & abs(markers$avg_log2FC) > lfc_cutoff] <- "Significant"
  
  # Select top N up and down regulated genes
  top_up <- markers %>%
    filter(avg_log2FC > 0, p_val_adj < pval_cutoff) %>%
    arrange(p_val_adj) %>%
    head(n = top_n)
  
  top_down <- markers %>%
    filter(avg_log2FC < 0, p_val_adj < pval_cutoff) %>%
    arrange(p_val_adj) %>%
    head(n = top_n)
  
  top_genes <- rbind(top_up, top_down)
  
  # Create plot
  p <- ggplot(markers, aes(x = avg_log2FC, y = log_pval, color = highlight, label = gene)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = "Log2 Fold Change",
      y = "-log10(Adjusted p-value)",
      color = "Significant",
      title = title
    )
  
  # Add labels
  if (nrow(top_genes) > 0) {
    if (label_all) {
      p <- p + geom_text_repel(
        data = top_genes,
        size = 3,
        box.padding = 0.5,
        max.overlaps = Inf
      )
    } else {
      p <- p + geom_text_repel(
        data = top_genes,
        size = 3,
        box.padding = 0.5,
        max.overlaps = 20
      )
    }
  }
  
  return(p)
}

# Function to create a heatmap of differentially expressed genes
create_heatmap <- function(seurat_obj,
                         markers,
                         group.by = "seurat_clusters",
                         top_n = 20,
                         scale = "row",
                         show_rownames = TRUE,
                         fontsize_row = 8) {
  
  # Get top markers
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = abs(avg_log2FC))
  
  # Get expression data
  expr <- GetAssayData(seurat_obj, slot = "scale.data")
  
  # Subset to marker genes
  expr <- expr[unique(top_markers$gene), ]
  
  # Create annotation
  annotation_col <- data.frame(
    Group = seurat_obj@meta.data[[group.by]]
  )
  rownames(annotation_col) <- colnames(seurat_obj)
  
  # Create color palette
  colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(255)
  
  # Create heatmap
  pheatmap(
    mat = expr,
    scale = scale,
    show_rownames = show_rownames,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    color = colors,
    fontsize_row = fontsize_row,
    clustering_method = "complete",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean"
  )
}

# Function to perform pseudobulk analysis
pseudobulk_analysis <- function(seurat_obj,
                              group.by = "seurat_clusters",
                              sample.by = "orig.ident",
                              design = ~ group,
                              test = "DESeq2") {
  
  # Aggregate counts
  counts <- AggregateExpression(
    seurat_obj,
    group.by = c(sample.by, group.by),
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA
  
  # Create metadata
  metadata <- data.frame(
    sample = colnames(counts),
    group = gsub("^[^_]*_", "", colnames(counts))
  )
  
  if (test == "DESeq2") {
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = round(counts),
      colData = metadata,
      design = design
    )
    
    # Filter lowly expressed genes
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds)
    
    return(list(
      results = as.data.frame(res),
      dds = dds
    ))
  } else if (test == "limma") {
    # Convert to voom object
    v <- voom(counts, design = model.matrix(design, metadata))
    
    # Fit linear model
    fit <- lmFit(v, design = model.matrix(design, metadata))
    fit <- eBayes(fit)
    
    # Get results
    res <- topTable(fit, number = Inf, sort.by = "none")
    
    return(list(
      results = res,
      voom = v,
      fit = fit
    ))
  }
}

# Main function for differential expression analysis
perform_differential_expression <- function(seurat_obj,
                                          group.by = "seurat_clusters",
                                          ident.1 = NULL,
                                          ident.2 = NULL,
                                          test.use = "wilcox",
                                          min.pct = 0.1,
                                          logfc.threshold = 0.25,
                                          only.pos = FALSE,
                                          top_n_genes = 10,
                                          output_dir = "results/differential_expression") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Find markers
  message("Finding differentially expressed genes...")
  markers <- find_markers(
    seurat_obj = seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    group.by = group.by,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos
  )
  
  # Save markers
  write.csv(
    markers,
    file.path(output_dir, "differential_genes.csv"),
    row.names = FALSE
  )
  
  # Create volcano plot
  message("Creating volcano plot...")
  volcano_plot <- create_volcano_plot(
    markers = markers,
    title = paste("Differential Expression:", ident.1, "vs", ifelse(is.null(ident.2), "Rest", ident.2)),
    top_n = top_n_genes
  )
  
  ggsave(
    filename = file.path(output_dir, "volcano_plot.png"),
    plot = volcano_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Create heatmap
  if (nrow(markers) > 0) {
    message("Creating heatmap...")
    png(
      filename = file.path(output_dir, "heatmap.png"),
      width = 10,
      height = 8,
      units = "in",
      res = 300
    )
    
    heatmap_plot <- create_heatmap(
      seurat_obj = seurat_obj,
      markers = markers,
      group.by = group.by,
      top_n = min(50, nrow(markers))
    )
    
    dev.off()
  } else {
    heatmap_plot <- NULL
  }
  
  # Return results
  return(list(
    markers = markers,
    volcano_plot = volcano_plot,
    heatmap_plot = heatmap_plot
  ))
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Example usage
  # library(Seurat)
  # 
  # # Load data
  # seurat_obj <- readRDS("results/clustering/seurat_clustered.rds")
  # 
  # # Run differential expression analysis
  # de_result <- perform_differential_expression(
  #   seurat_obj = seurat_obj,
  #   group.by = "seurat_clusters",
  #   ident.1 = "1",
  #   ident.2 = "2",
  #   test.use = "wilcox",
  #   output_dir = "results/differential_expression"
  # )
  # 
  # print("Differential expression analysis complete!")
}
