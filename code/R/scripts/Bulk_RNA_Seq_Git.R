# Load necessary libraries
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(pheatmap)
library(limma)
library(AnnotationDbi)
library(viridis)
library(grid) # for modifying pheatmap legend title
library(gtable)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(signatureSearch)

# -------------------------------
# Set base directory and data paths
# -------------------------------
base_dir <- "~/zachary-wilkes/scRNAseq_Publication/"

# Function to generate file paths
generate_file_paths <- function(base_dir) {
  list(
    "2D_1" = file.path(base_dir, "2D-01_S79_R1_001_cat", "2D-01_S79_R1_001_cat_trimmedReadsPerGene.out.tab"),
    "2D_2" = file.path(base_dir, "2D-02_S80_R1_001_cat", "2D-02_S80_R1_001_cat_trimmedReadsPerGene.out.tab"),
    "2D_3" = file.path(base_dir, "2D-03_S81_R1_001_cat", "2D-03_S81_R1_001_cat_trimmedReadsPerGene.out.tab"),
    "3D_1" = file.path(base_dir, "3D-01_S82_R1_001_cat", "3D-01_S82_R1_001_cat_trimmedReadsPerGene.out.tab"),
    "3D_2" = file.path(base_dir, "3D-02_S83_R1_001_cat", "3D-02_S83_R1_001_cat_trimmedReadsPerGene.out.tab"),
    "3D_3" = file.path(base_dir, "3D-03_S84_R1_001_cat", "3D-03_S84_R1_001_cat_trimmedReadsPerGene.out.tab"),
    "Fresh_1" = file.path(base_dir, "SRR9036929_Aligned", "SRR9036929_trimmedReadsPerGene.out.tab"),
    "Fresh_2" = file.path(base_dir, "SRR9036930_Aligned", "SRR9036930_trimmedReadsPerGene.out.tab"),
    "Fresh_3" = file.path(base_dir, "SRR9036931_Aligned", "SRR9036931_trimmedReadsPerGene.out.tab"),
    "Fresh_4" = file.path(base_dir, "SRR9036932_Aligned", "SRR9036932_trimmedReadsPerGene.out.tab")
  )
}

# Function to load and process data
load_and_process_data <- function(file_paths) {
  # Load files into a named list
  data_list <- lapply(file_paths, function(path) {
    read.table(path, sep = "\t", row.names = 1)
  })
  
  # Select relevant columns for each condition
  # 2D and 3D samples use V3, Fresh samples use V2
  sample_groups <- list(
    twoD_threeD = 1:6,
    fresh = 7:10
  )
  
  data_list[sample_groups$twoD_threeD] <- lapply(data_list[sample_groups$twoD_threeD], function(x) x[, "V3", drop = FALSE])
  data_list[sample_groups$fresh] <- lapply(data_list[sample_groups$fresh], function(x) x[, "V2", drop = FALSE])
  
  # Combine data into a single dataframe
  counts_df <- bind_cols(data_list)
  colnames(counts_df) <- names(data_list)
  
  # Remove unnecessary rows (first 4 rows are unwanted headers)
  counts_df <- counts_df[-c(1:4), ]
  
  # Ensure numeric format
  counts <- as.data.frame(lapply(counts_df, as.integer))
  rownames(counts) <- rownames(counts_df)
  colnames(counts) <- names(data_list)
  
  return(counts)
}

# Function to create metadata
create_metadata <- function(counts_df) {
  # Get sample counts for each condition
  num_2D <- sum(startsWith(colnames(counts_df), "2D"))
  num_3D <- sum(startsWith(colnames(counts_df), "3D"))
  num_Fresh <- sum(startsWith(colnames(counts_df), "Fresh"))
  
  # Create condition factor
  condition <- factor(c(rep("2D", num_2D), rep("3D", num_3D), rep("Fresh", num_Fresh)))
  
  # Create metadata dataframe
  metadata <- data.frame(Condition = condition)
  rownames(metadata) <- colnames(counts_df)
  
  return(metadata)
}

# Function to convert Ensembl IDs to gene symbols
convert_to_gene_symbols <- function(counts_df) {
  # Step 1: Remove version numbers from Ensembl IDs
  clean_ensembl_ids <- sub("\\..*$", "", rownames(counts_df))
  
  # Step 2: Map Ensembl IDs to gene symbols using org.Mm.eg.db
  gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys = clean_ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Step 3: Handle unmapped genes
  mapped <- !is.na(gene_symbols)
  counts_df <- counts_df[mapped, ]
  gene_symbols <- gene_symbols[mapped]
  
  # Step 4: Add gene symbols as a new column for aggregation
  counts_df <- counts_df %>%
    rownames_to_column(var = "Original_Ensembl_ID") %>%
    mutate(GeneSymbol = gene_symbols)
  
  # Step 5: Aggregate counts by GeneSymbol to handle duplicates
  counts_df <- counts_df %>%
    group_by(GeneSymbol) %>%
    summarise(across(starts_with("2D_") | starts_with("3D_") | starts_with("Fresh_"), sum, na.rm = TRUE)) %>%
    ungroup()
  
  # Step 6: Set GeneSymbol as row names
  counts_df <- counts_df %>%
    column_to_rownames(var = "GeneSymbol")
  
  # Verify that row names are now unique
  if(any(duplicated(rownames(counts_df)))) {
    stop("There are still duplicate gene symbols after aggregation.")
  }
  
  # Optional: Summary of mapping
  total_genes_initial <- length(clean_ensembl_ids)
  mapped_genes <- nrow(counts_df)
  unmapped_genes <- total_genes_initial - sum(mapped)
  
  cat("Total genes initially:", total_genes_initial, "\n")
  cat("Mapped genes:", mapped_genes, "\n")
  cat("Unmapped genes removed:", unmapped_genes, "\n")
  
  return(counts_df)
}

# Function to perform DESeq2 analysis
perform_deseq2_analysis <- function(counts_df, metadata) {
  # Create DESeq dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_df,
    colData = metadata,
    design = ~ Condition
  )
  
  # Filter low count genes
  keep <- rowSums(counts(dds) >= 60) >= 3
  dds <- dds[keep,]
  
  # Run DESeq
  dds <- DESeq(dds)
  
  return(dds)
}

# Function to extract and process DESeq2 results
extract_deseq2_results <- function(dds, condition1, condition2, output_path = NULL) {
  # Extract results for the comparison
  res <- results(dds, contrast = c("Condition", condition1, condition2))
  
  # Shrink log2 fold changes
  res <- lfcShrink(dds, contrast = c("Condition", condition1, condition2), res = res, type = "ashr")
  
  # Convert to data frame and add gene names as first column
  res_df <- as.data.frame(res)
  res_df <- add_column(res_df, rownames(res), .before = 1)
  colnames(res_df)[1] <- "gene"
  
  # Write to CSV if path is provided
  if (!is.null(output_path)) {
    write_csv(res_df, output_path)
  }
  
  return(res_df)
}

extract_deseq2_results_with_stat <- function(dds, condition1, condition2, output_path = NULL) {
  # Extract results for the comparison
  res <- results(dds, contrast = c("Condition", condition1, condition2))
  
  # Shrink log2 fold changes using ashr
  res <- lfcShrink(dds, contrast = c("Condition", condition1, condition2), res = res, type = "ashr")
  
  # Convert to data frame and add gene names as first column
  res_df <- as.data.frame(res)
  res_df <- add_column(res_df, rownames(res), .before = 1)
  colnames(res_df)[1] <- "gene"
  
  # Compute the test statistic: shrunken log2FoldChange divided by lfcSE
  res_df$stat <- res_df$log2FoldChange / res_df$lfcSE
  
  # Write to CSV if an output path is provided
  if (!is.null(output_path)) {
    write_csv(res_df, output_path)
  }
  
  return(res_df)
}

# Function to prepare gene list for GSEA
prepare_gene_list <- function(res) {
  if (!"stat" %in% colnames(res)) {
    if ("lfcSE" %in% colnames(res) && "log2FoldChange" %in% colnames(res)) {
      res$stat <- res$log2FoldChange / res$lfcSE
    } else {
      stop("DESeq2 result must include either stat or both log2FoldChange and lfcSE columns.")
    }
  }
  gene_list <- sort(na.omit(res$stat), decreasing = TRUE)
  names(gene_list) <- rownames(res)
  return(gene_list)
}

# Function to run GSEA GO for a specific ontology
run_gsea_go <- function(gene_list, ontology) {
  gseGO2(
    geneList = gene_list,
    ont = ontology,
    keyType = "SYMBOL",
    OrgDb = org.Mm.eg.db,
    minGSSize = 10,
    maxGSSize = 1000,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    pAdjustMethod = "BH",
    nPerm = 100000,
    nproc = 8
  )
}

# Function to simplify GO results
##simplify_go_results <- function(go_result) {
##  clusterProfiler::simplify(
##    go_result,
##  cutoff = 0.4,
##  by = "p.adjust",
##  select_fun = min,
##  measure = "Wang",
##  semData = NULL
##)
##}



# Function to create dotplot for GO results
create_go_dotplot <- function(go_result, title) {
  clusterProfiler::dotplot(go_result, showCategory = 20) + ggtitle(title)
}

# Function to generate heatmap of PCA distances
generate_pca_distance_heatmap <- function(sampleDistMatrix) {
  # Ensure row and column names are unique and correspond to samples
  if (anyDuplicated(rownames(sampleDistMatrix)) > 0) {
    stop("Non-unique row names detected in sampleDistMatrix. Please ensure each sample has a unique name.")
  }
  if (anyDuplicated(colnames(sampleDistMatrix)) > 0) {
    stop("Non-unique column names detected in sampleDistMatrix. Please ensure each sample has a unique name.")
  }
  
  # Define a suitable color palette for the heatmap
  heatmap_colors <- viridis(256)
  
  # Create the Heatmap object
  heatmap_obj <- Heatmap(sampleDistMatrix, 
                         cluster_rows = TRUE,    
                         cluster_columns = TRUE,  
                         show_row_names = TRUE,   
                         show_column_names = TRUE,
                         col = heatmap_colors,    
                         show_column_dend = TRUE,
                         show_row_dend = TRUE,
                         row_title = NULL,        
                         height = unit(5.5, "inch"),
                         heatmap_legend_param = list(
                           title = "PCA Distance",
                           title_gp = gpar(fontsize = 12, fontface = "bold"),
                           labels_gp = gpar(fontsize = 10), 
                           legend_width = unit(8, "cm"), 
                           legend_direction = "vertical",
                           just = "center"
                         )
  )
  
  # Draw the Heatmap and capture it as a grob
  heatmap_grob <- grid.grabExpr(draw(heatmap_obj, 
                                     heatmap_legend_side = "right",
                                     annotation_legend_side = "bottom",
                                     padding = unit(c(10, 10, 10, 10), "mm")))
  
  return(heatmap_grob)
}

# Function to create volcano plot
create_volcano_plot <- function(res, title) {
  # Convert DESeq2 results to a data frame and remove rows with NA in padj or log2FoldChange
  res_df <- res
  res_df <- res_df[complete.cases(res_df$log2FoldChange, res_df$padj), ]
  
  # Check if any genes remain after filtering
  if (nrow(res_df) == 0) {
    stop("No genes remain after filtering out NA values.")
  }
  
  # Ensure rownames exist and are unique
  if (is.null(rownames(res_df)) || anyDuplicated(rownames(res_df)) > 0) {
    stop("The results data frame must have unique row names.")
  }
  
  # Default all genes to grey (non-significant)
  res_df$col <- "grey"
  
  # Define significance criteria
  pcut <- 1e-05
  fccut <- 1.0
  
  pos_genes <- which(res_df$padj < pcut & res_df$log2FoldChange >= fccut)
  neg_genes <- which(res_df$padj < pcut & res_df$log2FoldChange <= -fccut)
  pval_only_genes <- which(res_df$padj < pcut & abs(res_df$log2FoldChange) < fccut)
  
  # Assign colors
  res_df$col[pval_only_genes] <- "blue"
  res_df$col[neg_genes] <- "red"
  res_df$col[pos_genes] <- "green"
  
  custom_col <- res_df$col
  
  names(custom_col)[custom_col == 'grey'] <- 'NS'
  names(custom_col)[custom_col == 'blue'] <- 'P-Value Significant'
  names(custom_col)[custom_col == 'red'] <- 'Significant and Negative'
  names(custom_col)[custom_col == 'green'] <- 'Significant and Positive'
  
  # Get top genes
  top_up <- res_df %>%
    mutate(gene = rownames(res_df)) %>%
    arrange(desc(log2FoldChange)) %>%
    head(20) %>%
    pull(gene)
  
  top_down <- res_df %>%
    mutate(gene = rownames(res_df)) %>%
    arrange(log2FoldChange) %>%
    head(20) %>%
    pull(gene)
  
  top_genes <- c(top_up, top_down)
  
  # Create volcano plot
  p <- EnhancedVolcano(
    res_df,
    lab = rownames(res_df),
    selectLab = top_genes,
    x = 'log2FoldChange',
    y = 'padj',
    title = title,
    arrowheads = FALSE,
    pCutoff = pcut,
    FCcutoff = fccut,
    pointSize = 1,
    subtitle = NULL,
    xlim = c(-20, 20),
    ylim = c(0, 326),
    labSize = 3,
    colCustom = custom_col,
    colAlpha = 0.6,
    drawConnectors = TRUE,
    max.overlaps = Inf,
    widthConnectors = 0.5,
    colConnectors = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = "full",
    borderWidth = 0.5,
    borderColour = "black",
    caption = NULL) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  return(p)
}

# Function to save MA plot
save_ma_plot <- function(res, filename, title) {
  png(filename, width = 6, height = 6, units = "in", res = 300)
  plotMA(res, main = title)
  dev.off()
}

# -------------------------------
# MAIN WORKFLOW
# -------------------------------

# 1. Generate file paths and load data
file_paths <- generate_file_paths(base_dir)
counts_df <- load_and_process_data(file_paths)

# 2. Create metadata
metadata <- create_metadata(counts_df)

# 3. Convert Ensembl IDs to gene symbols
counts_df <- convert_to_gene_symbols(counts_df)

# 4. Perform DESeq2 analysis
dds <- perform_deseq2_analysis(counts_df, metadata)

# 5. Extract DESeq2 results for all comparisons
comparisons <- list(
  "2D_vs_Fresh" = list(condition1 = "2D", condition2 = "Fresh", 
                       output_path = "~/zachary-wilkes/scRNAseq_Publication/res_2D_vs_Fresh_df.csv"),
  "3D_vs_Fresh" = list(condition1 = "3D", condition2 = "Fresh", 
                       output_path = "~/zachary-wilkes/scRNAseq_Publication/res_3D_vs_Fresh_df.csv"),
  "2D_vs_3D" = list(condition1 = "2D", condition2 = "3D", 
                    output_path = "~/zachary-wilkes/scRNAseq_Publication/res_2D_vs_3D_df.csv")
)

results_list <- lapply(names(comparisons), function(comp_name) {
  comp <- comparisons[[comp_name]]
  result <- extract_deseq2_results(dds, comp$condition1, comp$condition2, comp$output_path)
  return(result)
})
names(results_list) <- names(comparisons)

# 6. Data transformation for visualization
vsd <- DESeq2::vst(dds, blind = FALSE)

# 7. PCA Plot
pca_plot <- plotPCA(vsd, intgroup = "Condition") + ggtitle("PCA Plot")
print(pca_plot)

# 8. Sample-to-sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- colnames(vsd)
rownames(sampleDistMatrix) <- colnames(vsd)

annotation_col <- data.frame(Condition = vsd$Condition)
rownames(annotation_col) <- colnames(vsd)

heatmap_Bulk_grob <- generate_pca_distance_heatmap(sampleDistMatrix)
heatmap_Bulk_title <- ggdraw() + 
  draw_label("Heatmap of PCA Distance Matrices", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_bw() +
  theme(plot.margin = margin(-19, 35, 0, 0))

heatmap_Bulk_plot <- plot_grid(heatmap_Bulk_title, heatmap_Bulk_grob, ncol = 1, rel_heights = c(0.1, 1))
ggsave("heatmap_Bulk_plot.png", heatmap_Bulk_plot, width = 8, height = 8, dpi = 300)

# 9. Generate volcano plots
volcano_plots <- lapply(names(results_list), function(comp_name) {
  create_volcano_plot(results_list[[comp_name]], comp_name)
})
names(volcano_plots) <- names(results_list)

# Save individual volcano plots
ggsave("volcano_2D_vs_Fresh.png", volcano_plots[["2D_vs_Fresh"]], width = 8, height = 8, dpi = 300)

# Print volcano plots
for(plot in volcano_plots) {
  print(plot)
}

# Combined volcano plot
combined_volcano <- plot_grid(
  volcano_plots[["2D_vs_Fresh"]],
  volcano_plots[["3D_vs_Fresh"]],
  volcano_plots[["2D_vs_3D"]],
  labels = c("A", "B", "C"),
  ncol = 3,
  align = "hv",
  axis = "l"
)


ggsave("combined_volcano_figure.png", combined_volcano, width = 24, height = 8, dpi = 300)

# 10. Prepare gene lists for GSEA
gene_lists <- lapply(names(results_list), function(comp_name) {
  prepare_gene_list(results_list[[comp_name]])
})
names(gene_lists) <- names(results_list)

# 11. Run GSEA GO for all comparisons and ontologies
ontologies <- c("BP", "CC", "MF")
go_results <- list()

for(comp_name in names(gene_lists)) {
  go_results[[comp_name]] <- list()
  for(ont in ontologies) {
    go_result <- run_gsea_go(gene_lists[[comp_name]], ont)
    #go_result <- simplify_go_results(go_result)
    go_results[[comp_name]][[ont]] <- result(go_result)
  }
}

# 12. Create dotplots for GO results
dotplots <- list()

for(comp_name in names(go_results)) {
  dotplots[[comp_name]] <- list()
  for(ont in ontologies) {
    title <- paste0("GO ", 
                    ifelse(ont == "BP", "Biological Process", 
                           ifelse(ont == "CC", "Cellular Component", "Molecular Function")),
                    ": ", comp_name)
    dotplots[[comp_name]][[ont]] <- create_go_dotplot(go_results[[comp_name]][[ont]], title)
  }
}

# 13. Combine dotplots for each comparison
combined_dotplots <- list()

for(comp_name in names(dotplots)) {
  combined_dotplots[[comp_name]] <- plot_grid(
    dotplots[[comp_name]][["BP"]],
    dotplots[[comp_name]][["CC"]],
    dotplots[[comp_name]][["MF"]],
    labels = c("A", "B", "C"),
    ncol = 3,
    align = "hv",
    axis = "l"
  )
  
  # Save combined dotplots
  ggsave(paste0("combined_GO_ontologies_", comp_name, ".png"), 
         combined_dotplots[[comp_name]], 
         width = 24, height = 10, dpi = 300)
  
  print(combined_dotplots[[comp_name]])
}

# 14. Combine all comparison dotplots
combined_GO_ontologies <- plot_grid(
  combined_dotplots[["2D_vs_Fresh"]],
  combined_dotplots[["3D_vs_Fresh"]],
  combined_dotplots[["2D_vs_3D"]],
  nrow = 3,
  ncol = 1,
  align = "hv",
  axis = "l"
)

ggsave("combined_GO_ontologies.png", combined_GO_ontologies, width = 20, height = 28, dpi = 300)

# 15. Create MA plots
for(comp_name in names(results_list)) {
  save_ma_plot(
    results_list[[comp_name]], 
    paste0("MA_", comp_name, ".png"),
    paste0("MA Plot: ", comp_name)
  )
}

# 16. Hierarchical clustering dendrogram of samples
hc_samples <- hclust(dist(t(assay(dds))))
png("sample_clustering_dendrogram.png", width = 6, height = 6, units = "in", res = 300)
plot(hc_samples, main = "Hierarchical Clustering of Samples", xlab = "", sub = "")
dev.off()

# 17. Dispersion estimates plot
png("dispersion_estimates.png", width = 6, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()
