# ============================================
# Integrated Analysis: Mouse RNA-seq & Proteomics
# ============================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# --------------------------------------------
# 1. Data Loading
# --------------------------------------------
# Replace the file paths with your actual data files.
# Each file is assumed to include at least:
#   - A common gene identifier: 'gene_id'
#   - A log2 fold change column: 'log2FC'
#   - A p-value column: 'pvalue'

bulk_results <- read.csv("bulk_results.csv", stringsAsFactors = FALSE)
sc_results   <- read.csv("scRNAseq_results.csv", stringsAsFactors = FALSE)
prot_results <- read.csv("proteomics_results.csv", stringsAsFactors = FALSE)

# --------------------------------------------
# 2. Merge Datasets for Correlation Analysis
# --------------------------------------------
# Merge bulk RNA-seq and proteomics data by 'gene_id'
bulk_prot <- inner_join(bulk_results, prot_results, by = "gene_id", 
                        suffix = c("_bulk", "_prot"))

# Merge single cell RNA-seq and proteomics data by 'gene_id'
sc_prot <- inner_join(sc_results, prot_results, by = "gene_id", 
                      suffix = c("_sc", "_prot"))

# --------------------------------------------
# 3. Correlation Analyses
# --------------------------------------------
# Pearson correlation: Bulk RNA-seq vs. Proteomics
bulk_cor <- cor.test(bulk_prot$log2FC_bulk, bulk_prot$log2FC_prot, method = "pearson")
cat("Bulk RNA-seq vs. Proteomics correlation:\n")
print(bulk_cor)

# Pearson correlation: Single Cell RNA-seq vs. Proteomics
sc_cor <- cor.test(sc_prot$log2FC_sc, sc_prot$log2FC_prot, method = "pearson")
cat("\nSingle Cell RNA-seq vs. Proteomics correlation:\n")
print(sc_cor)

# --------------------------------------------
# 4. Visualization: Scatter Plots with Regression Lines
# --------------------------------------------
# Scatter plot: Bulk RNA-seq vs. Proteomics
p_bulk <- ggplot(bulk_prot, aes(x = log2FC_bulk, y = log2FC_prot)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Bulk RNA-seq vs. Proteomics",
       x = "Bulk RNA log2 Fold Change",
       y = "Proteomics log2 Fold Change") +
  theme_minimal()

# Scatter plot: Single Cell RNA-seq vs. Proteomics
p_sc <- ggplot(sc_prot, aes(x = log2FC_sc, y = log2FC_prot)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "scRNA-seq vs. Proteomics",
       x = "scRNA-seq log2 Fold Change",
       y = "Proteomics log2 Fold Change") +
  theme_minimal()

# Print plots to the R graphics device
print(p_bulk)
print(p_sc)

# --------------------------------------------
# 5. Overlap Analysis: Venn Diagram
# --------------------------------------------
# Identify significant genes based on a chosen threshold (e.g., pvalue < 0.05)
sig_sc   <- unique(sc_results[sc_results$pvalue < 0.05, "gene_id"])
sig_prot <- unique(prot_results[prot_results$pvalue < 0.05, "gene_id"])

# Create a Venn diagram to show the overlap between scRNA-seq and proteomics significant genes
venn.diagram(
  x = list("scRNA-seq" = sig_sc, "Proteomics" = sig_prot),
  filename = "venn_diagram.png",
  fill = c("blue", "green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap between scRNA-seq and Proteomics"
)
cat("\nVenn diagram saved as 'venn_diagram.png'\n")

# --------------------------------------------
# 6. Pathway Enrichment Analysis
# --------------------------------------------
# For mouse data, we convert gene symbols to Entrez IDs using org.Mm.eg.db.
# Here we perform enrichment analysis using KEGG pathways.

# Extract significant genes (using the same threshold as above)
sig_sc_genes   <- unique(sc_results[sc_results$pvalue < 0.05, "gene_id"])
sig_prot_genes <- unique(prot_results[prot_results$pvalue < 0.05, "gene_id"])

# Convert gene symbols to Entrez IDs
sc_entrez   <- bitr(sig_sc_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
prot_entrez <- bitr(sig_prot_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform KEGG enrichment analysis for scRNA-seq significant genes
enrich_sc <- enrichKEGG(gene         = sc_entrez$ENTREZID,
                        organism     = 'mmu',
                        pvalueCutoff = 0.05)

# Perform KEGG enrichment analysis for Proteomics significant genes
enrich_prot <- enrichKEGG(gene         = prot_entrez$ENTREZID,
                          organism     = 'mmu',
                          pvalueCutoff = 0.05)

cat("\nTop enriched KEGG pathways for scRNA-seq data:\n")
print(head(enrich_sc))

cat("\nTop enriched KEGG pathways for Proteomics data:\n")
print(head(enrich_prot))

# Visualize the top 10 enriched pathways with barplots
barplot(enrich_sc, showCategory = 10, title = "KEGG Enrichment: scRNA-seq")
barplot(enrich_prot, showCategory = 10, title = "KEGG Enrichment: Proteomics")

# ============================================
# End of Analysis Script
# ============================================
