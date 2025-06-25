if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
options(repos = BiocManager::repositories())

install_if_missing  <- function(pkg) if (!requireNamespace(pkg, quietly = TRUE))
  install.packages(pkg, dependencies = TRUE)

install_if_missing_bioc <- function(pkg) if (!requireNamespace(pkg, quietly = TRUE))
  BiocManager::install(pkg, ask = FALSE)


install_if_missing_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE))
    devtools::install_github(repo, dependencies = TRUE)  # now pulls BioC deps
}

remotes::install_github("r-spatial/s2@master", dependencies = TRUE)
## ---- 1. Ensure Bioc + devtools ---------------------------------------------
install_if_missing("BiocManager")
install_if_missing("devtools")       # brings in remotes, pkgbuild, etc.

## ---- 2. Define package sets -----------------------------------------------
bioc_pkgs   <- c("AnnotationHub","biomaRt","clusterProfiler","decoupleR",
                 "dorothea","enrichplot","ensembldb","KEGGREST","MAST",
                 "OmnipathR","scater","viper","zinbwave", "sparseMatrixStats", "progeny",
                 "singscore", "GSVA", "GSEABase")

cran_pkgs   <- c("bigmemory","cowplot","data.table","dplyr","DT","flexmix",
                 "flextable","foreach","future","ggplot2","ggpubr","igraph",
                 "KernSmooth","kableExtra","knitr","magrittr","msigdbr",
                 "patchwork","pheatmap","progeny","RCurl","RColorBrewer",
                 "ROCR","SAVER","sctransform","Seurat","stringi","stringr",
                 "tibble","tidyr","uwot","writexl","xlsx","checkmate",
                 "clustermole")

gh_pkgs <- c("SeuratWrappers"="satijalab/seurat-wrappers",
             "Stringendo"       = "vertesy/Stringendo",                  # string helpers :contentReference[oaicite:0]{index=0}
             "CodeAndRoll2"     = "vertesy/CodeAndRoll2",                # 130+ utility funcs :contentReference[oaicite:1]{index=1}
             "ReadWriter"       = "vertesy/ReadWriter",                  # I/O convenience :contentReference[oaicite:2]{index=2}
             "MarkdownHelpers"  = "vertesy/MarkdownHelpers",             # Rmd helpers :contentReference[oaicite:3]{index=3}
             "MarkdownReports"  = "vertesy/MarkdownReports",             # auto-reporting :contentReference[oaicite:4]{index=4}
             "ggExpress"        = "vertesy/ggExpress",                   # fast plotting :contentReference[oaicite:5]{index=5}
             "DatabaseLinke.R"  = "vertesy/DatabaseLinke.R",     
             "Seurat.utils"  ="vertesy/Seurat.utils",
             "SeuratData"    ="satijalab/seurat-data",
             "DoubletFinder" ="chris-mcginnis-ucsf/DoubletFinder",
             "Lamian"        ="Winnie09/Lamian",
             "Libra"         ="cran/Libra",
             "PISCES"        ="califano-lab/PISCES",
             "SCPA"          ="jackbibby1/SCPA",
             "colorblindr"   ="clauswilke/colorblindr",
             "multicross" = "cran/multicross")

## ---- 3. Install in priority order -----------------------------------------
lapply(bioc_pkgs,   install_if_missing_bioc)
lapply(cran_pkgs,   install_if_missing)
mapply(install_if_missing_github, names(gh_pkgs), gh_pkgs, USE.NAMES = FALSE)
remotes::install_version("qpcR", version = "1.4-1",
                         repos = "https://cran.r-project.org")


## ---- 4. Load everything (optional) ----------------------------------------
invisible(lapply(c(bioc_pkgs, cran_pkgs, names(gh_pkgs)),
                 function(pkg) suppressPackageStartupMessages(library(pkg,
                                                                      character.only = TRUE))))





setwd("~/scRNAseq_Publication/")

if(require("flexiblas", quietly = TRUE)) {
  flexiblas_load_backend("OPENBLAS-SERIAL")
  flexiblas_switch(n = 2)
  flexiblas_set_num_threads(1)
  message("Configured FlexiBLAS for parallel operations")
} else {
  message("FlexiBLAS not available. For best performance, consider installing it.")
}
#####
GFP_Ag_2D_LP.data <- Read10X(data.dir = "Data/aTomei_zWilkes_scRNAseq_12162021/multi_sample/GFP-Ag_2D_LP_sample_feature_bc_matrix")
GFP_Ag_3D_LP.data <- Read10X(data.dir = "Data/aTomei_zWilkes_scRNAseq_12162021/multi_sample/GFP-Ag_3D_LP_sample_feature_bc_matrix")
Fresh_FRCs.data <- Read10X(data.dir = "Data/aTomei_zWilkes_scRNAseq_12162021/single_sample/10X1028GeXSI_12162021_filtered_feature_bc_matrix")
#####
GFP_Ag_2D_LP <- CreateSeuratObject(counts = GFP_Ag_2D_LP.data$`Gene Expression`, project = "GFPAg2DLP", min.cells = 3, min.features = 200)
GFP_Ag_3D_LP <- CreateSeuratObject(counts = GFP_Ag_3D_LP.data$`Gene Expression`, project = "GFPAg3DLP", min.cells = 3, min.features = 200)
Fresh_FRCs <- CreateSeuratObject(counts = Fresh_FRCs.data, project = "FreshFRCs", min.cells = 3, min.features = 200)
#####
#FRC.Combined.FreshVCulture <- merge(GFP_Ag_2D_LP, y = c(GFP_Ag_3D_LP, Fresh_FRCs)
#                                    , add.cell.ids = c("GFPAg2DLP", "GFPAg3DLP", "FreshlyIsolated"), project = "FRC.combined.FreshVCulture")
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s.genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m.genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


GFP_Ag_2D_LP$log10GenesPerUMI <- log10(GFP_Ag_2D_LP$nFeature_RNA) / log10(GFP_Ag_2D_LP$nCount_RNA)
GFP_Ag_3D_LP$log10GenesPerUMI <- log10(GFP_Ag_3D_LP$nFeature_RNA) / log10(GFP_Ag_3D_LP$nCount_RNA)
Fresh_FRCs$log10GenesPerUMI <- log10(Fresh_FRCs$nFeature_RNA) / log10(Fresh_FRCs$nCount_RNA)


temp_list <- list(GFP_Ag_2D_LP, GFP_Ag_3D_LP, Fresh_FRCs)
temp_list <- setNames(temp_list, c("GFP_Ag_2D_LP", "GFP_Ag_3D_LP", "Fresh_FRCs"))
empty_list <- list()

for (i in temp_list) {
  i[["percent.mt"]] <- PercentageFeatureSet(i, pattern = "^mt-")
  empty_list <- append(empty_list, i)
}

temp_list <- empty_list
empty_list <- list()
#####

violin_plot = function(FRC_input){
  VlnPlot(FRC_input, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), stack = T) + 
    geom_point(position = position_jitter(seed = 1, width = 0.15), size = 0.2) 
}

scatterplot_ncount_nfeature = function(FRC_input){
  FeatureScatter(FRC_input, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.55)
}

scatterplot_ncount_percent.mt = function(FRC_input){
  FeatureScatter(FRC_input, feature1 = "nCount_RNA" , feature2 = "percent.mt", pt.size = 0.55)
}

#####

violin_plots <- lapply(temp_list, violin_plot)
scatter_ncount_nfeature <- lapply(temp_list, scatterplot_ncount_nfeature)
scatter_ncount_percent.mt <- lapply(temp_list, scatterplot_ncount_percent.mt)

violin_plots[1]
violin_plots[2]
violin_plots[3]

scatter_ncount_nfeature[1]
scatter_ncount_nfeature[2]
scatter_ncount_nfeature[3]

scatter_ncount_percent.mt[1]
scatter_ncount_percent.mt[2]
scatter_ncount_percent.mt[3]

#####

for (i in temp_list) {
  i <- RunMiQC(i, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.95)
  empty_list <- append(empty_list, i)
}

temp_list <- empty_list
empty_list <- list()

PlotMiQC(temp_list[[1]], color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
PlotMiQC(temp_list[[1]], color.by = "miQC.keep")
PlotMiQC(temp_list[[2]], color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
PlotMiQC(temp_list[[2]], color.by = "miQC.keep")
PlotMiQC(temp_list[[3]], color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
PlotMiQC(temp_list[[3]], color.by = "miQC.keep")

for (i in temp_list) {
  i <- subset(i, miQC.keep == "keep")
  empty_list <- append(empty_list, i)
}

temp_list <- empty_list
empty_list <- list()

violin_plots_2 <- lapply(temp_list, violin_plot)
scatter_ncount_nfeature_2 <- lapply(temp_list, scatterplot_ncount_nfeature)
scatter_ncount_percent.mt_2 <- lapply(temp_list, scatterplot_ncount_percent.mt)


violin_plots_2[1]
violin_plots_2[2]
violin_plots_2[3]

scatter_ncount_nfeature_2[1]
scatter_ncount_nfeature_2[2]
scatter_ncount_nfeature_2[3]

scatter_ncount_percent.mt_2[1]
scatter_ncount_percent.mt_2[2]
scatter_ncount_percent.mt_2[3]

#
for (i in temp_list) {
  i <- subset(i,  (nCount_RNA >= 500) & 
                (nFeature_RNA >= 250) )
  empty_list <- append(empty_list, i)
}

temp_list <- empty_list
empty_list <- list()

violin_plots_3 <- lapply(temp_list, violin_plot)
scatter_ncount_nfeature_3 <- lapply(temp_list, scatterplot_ncount_nfeature)
scatter_ncount_percent.mt_3 <- lapply(temp_list, scatterplot_ncount_percent.mt)


violin_plots_3[1]
violin_plots_3[2]
violin_plots_3[3]

scatter_ncount_nfeature_3[1]
scatter_ncount_nfeature_3[2]
scatter_ncount_nfeature_3[3]

scatter_ncount_percent.mt_3[1]
scatter_ncount_percent.mt_3[2]
scatter_ncount_percent.mt_3[3]

sctransform_fun = function(FRC_input){
  SCTransform(FRC_input, verbose = T)
}

temp_list <- lapply(temp_list, sctransform_fun)



temp_list <- setNames(temp_list, c("GFP_Ag_2D_LP", "GFP_Ag_3D_LP", "Fresh_FRCs"))
empty_list <- list()


PCA_fun = function(FRC_input){
  RunPCA(FRC_input, verbose = T, assay = "SCT")
}

temp_list <- lapply(temp_list, PCA_fun)

ElbowPlot(temp_list[[1]])
ElbowPlot(temp_list[[2]])
ElbowPlot(temp_list[[3]])

pc_list <- list()
for (i in 1:3) {
  stdv <- temp_list[[i]]@reductions$pca@stdev
  sum.stdv <- sum(temp_list[[i]]@reductions$pca@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  pc_list[i] <- min(co1, co2)
}

GFP_Ag_2D_LP <- temp_list[[1]] 
GFP_Ag_3D_LP <- temp_list[[2]]
Fresh_FRCs <- temp_list[[3]]


GFP_Ag_2D_LP <- RunUMAP(GFP_Ag_2D_LP, dims = 1:50, verbose = T, min.dist = 0.2, n.epochs = 5000, spread = 1, n.neighbors = 50, metric = "correlation", assay = "SCT")
GFP_Ag_2D_LP <- FindNeighbors(GFP_Ag_2D_LP, dims = 1:2, verbose = T, reduction = "umap", annoy.metric = "euclidean", n.trees = 3)
GFP_Ag_2D_LP <- FindClusters(GFP_Ag_2D_LP, verbose = T, resolution = 0.5)


GFP_Ag_3D_LP <- RunUMAP(GFP_Ag_3D_LP, dims = 1:50, verbose = T, min.dist = 0.2, n.epochs = 5000, spread = 1, n.neighbors = 50, metric = "correlation", assay = "SCT")
GFP_Ag_3D_LP <- FindNeighbors(GFP_Ag_3D_LP, dims = 1:2, verbose = T, reduction = "umap", annoy.metric = "euclidean", n.trees = 3)
GFP_Ag_3D_LP <- FindClusters(GFP_Ag_3D_LP, verbose = T, resolution = 0.5)


Fresh_FRCs <- RunUMAP(Fresh_FRCs, dims = 1:50, verbose = T, min.dist = 0.2, n.epochs = 5000, spread = 1, n.neighbors = 50, metric = "correlation", assay = "SCT")
Fresh_FRCs <- FindNeighbors(Fresh_FRCs, dims = 1:2, verbose = T, reduction = "umap", annoy.metric = "euclidean", n.trees = 3)
Fresh_FRCs <- FindClusters(Fresh_FRCs, verbose = T, resolution = 0.5)


DimPlot(GFP_Ag_2D_LP)
DimPlot(GFP_Ag_3D_LP)
DimPlot(Fresh_FRCs)


############################################################################################################################################################


process_and_find_doublets <- function(seurat_obj, pN = 0.10) {
  # Parameter sweep and optimization
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:50, sct = TRUE, num.cores = 16)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]  
  
  
  # Homotypic proportion calculation
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(seurat_obj@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # Doublet identification and filtering
  doublets <- doubletFinder(seurat_obj, PCs = 1:50, pN = pN, pK = optimal.pk, 
                            nExp = nExp.poi.adj, sct = TRUE)
  
  # Set identity based on doublet classifications
  doublets <- SetIdent(doublets, value = doublets@meta.data[[grep("^DF\\.classifications", colnames(doublets@meta.data), value = TRUE)]])
  
  # Visualize and subset singlets
  print(DimPlot(doublets))
  doublets <- subset(doublets, idents = "Singlet")
  print(DimPlot(doublets))
  
  return(doublets)
}

# Example Usage (assuming pc_list is defined elsewhere)
doublets_GFP_Ag_2D_LP <- process_and_find_doublets(GFP_Ag_2D_LP)
doublets_GFP_Ag_3D_LP <- process_and_find_doublets(GFP_Ag_3D_LP)
doublets_Fresh_FRCs    <- process_and_find_doublets(Fresh_FRCs)
#############################################################################################################################################################
doublets_GFP_Ag_2D_LP <- SetIdent(doublets_GFP_Ag_2D_LP, value = doublets_GFP_Ag_2D_LP$orig.ident)
doublets_GFP_Ag_3D_LP <- SetIdent(doublets_GFP_Ag_3D_LP, value = doublets_GFP_Ag_3D_LP$orig.ident)
doublets_Fresh_FRCs <- SetIdent(doublets_Fresh_FRCs, value = doublets_Fresh_FRCs$orig.ident)

GFP_Ag_2D_LP <- doublets_GFP_Ag_2D_LP
GFP_Ag_3D_LP <- doublets_GFP_Ag_3D_LP
Fresh_FRCs <- doublets_Fresh_FRCs

#########################################################################################################################################

cl <- parallel::makeCluster(8, type = "FORK")  # using 4 cores

doParallel::registerDoParallel(cl)
# Total number of predicted genes
total_genes <- length(GFP_Ag_2D_LP@assays$RNA@features)

# Compute 5 breakpoints that split 1:total_genes into 4 groups
breaks <- round(seq(1, total_genes + 1, length.out = 5))

# Create a list of 4 intervals (each interval is a vector of gene indices)
pred_ranges <- lapply(1:4, function(i) breaks[i]:(breaks[i + 1] - 1))
# For example, this will generate:
# pred_ranges[[1]] will be 1:3856, 
# pred_ranges[[2]] will be 3857:7712, etc.

# Now run the saver function for each gene subset
saver_list <- lapply(pred_ranges, function(range) {
  saver(
    GFP_Ag_2D_LP@assays$RNA$counts,
    pred.genes = range,
    pred.genes.only = TRUE,
    do.fast = TRUE
  )
})

# Optionally, name each element of the result for clarity:
names(saver_list) <- paste0("saver", 1:4, "_GFP_Ag_2D_LP")
saver.all_GFP_Ag_2D_LP <- combine.saver(saver_list)
GFP_Ag_2D_LP_SAVER <- CreateSeuratObject(counts = saver.all_GFP_Ag_2D_LP$estimate, assay = "SAVER")

GFP_Ag_2D_LP_SAVER[["percent.mt"]] <- PercentageFeatureSet(GFP_Ag_2D_LP_SAVER, pattern = "^mt-")

VlnPlot(GFP_Ag_2D_LP_SAVER, features = c("nFeature_SAVER", "nCount_SAVER", "percent.mt"), stack = T) + 
  geom_point(position = position_jitter(seed = 1, width = 0.15), size = 0.2) 

GFP_Ag_2D_LP[["SAVER"]] <- GFP_Ag_2D_LP_SAVER@assays$SAVER

FeatureScatter(GFP_Ag_2D_LP_SAVER, feature1 = "nCount_SAVER", feature2 = "nFeature_SAVER", pt.size = 0.55)



FeatureScatter(GFP_Ag_2D_LP_SAVER, feature1 = "nCount_SAVER" , feature2 = "percent.mt", pt.size = 0.55)





parallel::stopCluster(cl)
cl <- parallel::makeCluster(8, type = "FORK")  # using 4 cores
doParallel::registerDoParallel(cl)

# Total number of predicted genes
total_genes <- length(GFP_Ag_3D_LP@assays$RNA@features)

# Compute 5 breakpoints that split 1:total_genes into 4 groups
breaks <- round(seq(1, total_genes + 1, length.out = 5))

# Create a list of 4 intervals (each interval is a vector of gene indices)
pred_ranges <- lapply(1:4, function(i) breaks[i]:(breaks[i + 1] - 1))
# For example, this will generate:
# pred_ranges[[1]] will be 1:3856, 
# pred_ranges[[2]] will be 3857:7712, etc.

# Now run the saver function for each gene subset
saver_list <- lapply(pred_ranges, function(range) {
  saver(
    GFP_Ag_3D_LP@assays$RNA$counts,
    pred.genes = range,
    pred.genes.only = TRUE,
    do.fast = TRUE
  )
})

# Optionally, name each element of the result for clarity:
names(saver_list) <- paste0("saver", 1:4, "_GFP_Ag_3D_LP")
saver.all_GFP_Ag_3D_LP <- combine.saver(saver_list)
GFP_Ag_3D_LP_SAVER <- CreateSeuratObject(counts = saver.all_GFP_Ag_3D_LP$estimate, assay = "SAVER")

GFP_Ag_3D_LP_SAVER[["percent.mt"]] <- PercentageFeatureSet(GFP_Ag_3D_LP_SAVER, pattern = "^mt-")

GFP_Ag_3D_LP[["SAVER"]] <- GFP_Ag_3D_LP_SAVER@assays$SAVER

parallel::stopCluster(cl)
cl <- parallel::makeCluster(8, type = "FORK")  # using 4 cores
doParallel::registerDoParallel(cl)

# Total number of predicted genes
total_genes <- length(Fresh_FRCs@assays$RNA@features)

# Compute 5 breakpoints that split 1:total_genes into 4 groups
breaks <- round(seq(1, total_genes + 1, length.out = 5))

# Create a list of 4 intervals (each interval is a vector of gene indices)
pred_ranges <- lapply(1:4, function(i) breaks[i]:(breaks[i + 1] - 1))
# For example, this will generate:
# pred_ranges[[1]] will be 1:3856, 
# pred_ranges[[2]] will be 3857:7712, etc.

# Now run the saver function for each gene subset
saver_list <- lapply(pred_ranges, function(range) {
  saver(
    Fresh_FRCs@assays$RNA$counts,
    pred.genes = range,
    pred.genes.only = TRUE,
    do.fast = TRUE
  )
})

# Optionally, name each element of the result for clarity:
names(saver_list) <- paste0("saver", 1:4, "_Fresh_FRCs")
saver.all_Fresh_FRCs <- combine.saver(saver_list)
Fresh_FRCs_SAVER <- CreateSeuratObject(counts = saver.all_Fresh_FRCs$estimate, assay = "SAVER")

Fresh_FRCs_SAVER[["percent.mt"]] <- PercentageFeatureSet(Fresh_FRCs_SAVER, pattern = "^mt-")

Fresh_FRCs[["SAVER"]] <- Fresh_FRCs_SAVER@assays$SAVER

parallel::stopCluster(cl)
cl <- parallel::makeCluster(8, type = "FORK")  # using 4 cores
doParallel::registerDoParallel(cl)
gc()
#########################################################################################################################################
SaveSeuratRds(GFP_Ag_2D_LP, file="/home/zxw284/Pub_Figures/GFP_Ag_2D_LP_Obj")
GFP_Ag_2D_LP<-LoadSeuratRds(file="/home/zxw284/Pub_Figures/GFP_Ag_2D_LP_Obj")
SaveSeuratRds(GFP_Ag_3D_LP, file="/home/zxw284/Pub_Figures/GFP_Ag_3D_LP_Obj")
GFP_Ag_3D_LP<-LoadSeuratRds(file="/home/zxw284/Pub_Figures/GFP_Ag_3D_LP_Obj")
SaveSeuratRds(Fresh_FRCs, file="/home/zxw284/Pub_Figures/Fresh_FRCs_Paper_Obj")
Fresh_FRCs<-LoadSeuratRds(file="/home/zxw284/Pub_Figures/Fresh_FRCs_Paper_Obj")

temp_list <- list(GFP_Ag_2D_LP, GFP_Ag_3D_LP, Fresh_FRCs)
temp_list <- setNames(temp_list, c("GFP_Ag_2D_LP", "GFP_Ag_3D_LP", "Fresh_FRCs"))
empty_list <- list()

plan(multicore, workers = 8)
options(future.globals.maxSize = 200 * 1024^3)  # Increase allowed size to 50GB
options(future.rng.onMisuse = "ignore")        # Handle RNG issues in parallel

# Define the function to set assay AND return the object
set_default_saver <- function(seurat_obj) {
  # Check if the SAVER assay exists before setting it
  
  DefaultAssay(object = seurat_obj) <- "SAVER"
  
  return(seurat_obj) # <<<--- Crucial: Return the modified object
}

# Apply the function using lapply
# This creates a NEW list with the modified objects
temp_list_updated <- lapply(temp_list, set_default_saver)

# Optional: If you want to overwrite your original list
temp_list <- temp_list_updated

for (i in temp_list) {
  i[["percent.mt"]] <- PercentageFeatureSet(i, pattern = "^mt-")
  empty_list <- append(empty_list, i)
}

temp_list <- empty_list
empty_list <- list()

normalize_fun = function(FRC_input){
  NormalizeData(FRC_input, assay = "SAVER")
}

temp_list <- lapply(temp_list, normalize_fun)

vst_fun = function(FRC_input){
  FindVariableFeatures(FRC_input, assay = "SAVER")
}

temp_list <- lapply(temp_list, vst_fun)

initial_scale_fun = function(FRC_input){
  ScaleData(FRC_input, features = rownames(FRC_input) ,vars.to.regress = c("percent.mt"), do.scale = T, do.center = T, assay = "SAVER")
}

temp_list <- lapply(temp_list, initial_scale_fun)

Cell_Cycle_Scoring_fun = function(FRC_input){
  CellCycleScoring(FRC_input, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
}

temp_list <- lapply(temp_list, Cell_Cycle_Scoring_fun)


final_scale_reg_fun = function(FRC_input){
  ScaleData(FRC_input, features = rownames(FRC_input) ,vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), do.scale = T, do.center = T, assay = "SAVER")
}

temp_list <- lapply(temp_list, final_scale_reg_fun)

# Verify for one object (replace 'GFP_Ag_2D_LP' with the actual list name)

vst_fun = function(FRC_input){
  FindVariableFeatures(FRC_input, assay = "SAVER")
}

temp_list <- lapply(temp_list, vst_fun)

DimPlot(temp_list[[1]])
DimPlot(temp_list[[2]])
DimPlot(temp_list[[3]])


temp_list <- setNames(temp_list, c("GFP_Ag_2D_LP", "GFP_Ag_3D_LP", "Fresh_FRCs"))
empty_list <- list()



PCA_fun = function(FRC_input){
  RunPCA(FRC_input, verbose = T, assay = "SAVER")
}

temp_list <- lapply(temp_list, PCA_fun)

ElbowPlot(temp_list[[1]])
ElbowPlot(temp_list[[2]])
ElbowPlot(temp_list[[3]])

pc_list <- list()
for (i in 1:3) {
  stdv <- temp_list[[i]]@reductions$pca@stdev
  sum.stdv <- sum(temp_list[[i]]@reductions$pca@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  pc_list[i] <- min(co1, co2)
}

GFP_Ag_2D_LP <- temp_list[[1]] 
GFP_Ag_3D_LP <- temp_list[[2]]
Fresh_FRCs <- temp_list[[3]]



GFP_Ag_2D_LP <- RunUMAP(
  GFP_Ag_2D_LP, 
  dims = 1:25,
  umap.method = "umap-learn",
  verbose = TRUE, 
  min.dist = 0.2, 
  n.epochs = 5000, 
  spread     = 1, 
  n.neighbors = 25, 
  metric = "correlation", 
  assay = "SAVER",
  slot = "scale.data")

GFP_Ag_2D_LP <- FindNeighbors(GFP_Ag_2D_LP, dims = 1:2,verbose = T, reduction = "umap", annoy.metric = "euclidean", n.trees = 1000, compute.SNN = T, prune.SNN = 0, k.param = 25, nn.method = "annoy", nn.eps = 0, features = NULL, assay = NULL)
GFP_Ag_2D_LP <- FindClusters(GFP_Ag_2D_LP, verbose = T, resolution = 0.2, modularity.fxn = 1, algorithm = 2, n.start = 1000, n.iter = 1000, random.seed = 123)
DimPlot(GFP_Ag_2D_LP)
GFP_Ag_2D_LP_Markers <- FindAllMarkers(GFP_Ag_2D_LP, assay = "SCT", slot = "counts", test.use = "MAST", latent.vars = c("S.Score", "G2M.Score", "percent.mt"), min.diff.pct = 0.2)
write.csv(GFP_Ag_2D_LP_Markers, file = "/home/zachary-wilkes/scRNAseq_Publication/GFP_Ag_2D_LP_Markers")

GFP_Ag_3D_LP <- RunUMAP(
  GFP_Ag_3D_LP,
  umap.method = "umap-learn",
  dims = 1:25, 
  verbose = TRUE, 
  min.dist = 0.2, 
  n.epochs = 5000, 
  spread     = 1, 
  n.neighbors = 25, 
  metric = "correlation", 
  assay = "SAVER",
  slot = "scale.data"
)

GFP_Ag_3D_LP <- FindNeighbors(GFP_Ag_3D_LP, dims = 1:2,verbose = T, reduction = "umap", annoy.metric = "euclidean", n.trees = 1000, compute.SNN = T, prune.SNN = 0, k.param = 25, nn.method = "annoy", nn.eps = 0, features = NULL, assay = NULL)
GFP_Ag_3D_LP <- FindClusters(GFP_Ag_3D_LP, verbose = T, resolution = 0.2, modularity.fxn = 1, algorithm = 2, n.start = 1000, n.iter = 1000, random.seed = 123)
DimPlot(GFP_Ag_3D_LP)
GFP_Ag_3D_LP_Markers <- FindAllMarkers(GFP_Ag_3D_LP, assay = "SCT", slot = "counts", test.use = "MAST", latent.vars = c("S.Score", "G2M.Score", "percent.mt"), min.diff.pct = 0.2)
write.csv(GFP_Ag_3D_LP_Markers, file = "/home/zachary-wilkes/scRNAseq_Publication/GFP_Ag_3D_LP_Markers")

ElbowPlot(Fresh_FRCs, ndims = 43)


Fresh_FRCs <- RunUMAP(
  Fresh_FRCs,
  umap.method = "umap-learn", 
  dims = 1:43, 
  verbose = TRUE, 
  min.dist = 0.2, 
  n.epochs = 5000, 
  spread     = 1, 
  n.neighbors = 50, 
  metric = "correlation", 
  assay = "SAVER",
  slot = "scale.data"
)
Fresh_FRCs <- FindNeighbors(Fresh_FRCs, dims = 1:2,verbose = T, reduction = "umap", annoy.metric = "euclidean", n.trees = 1000, compute.SNN = T, prune.SNN = 0, k.param = 40, nn.method = "annoy", nn.eps = 0, features = NULL, assay = NULL)
Fresh_FRCs <- FindClusters(Fresh_FRCs, verbose = T, resolution = 0.4, modularity.fxn = 1, algorithm = 2, n.start = 1000, n.iter = 1000, random.seed = 123)
DimPlot(Fresh_FRCs)
Fresh_FRCs_Markers <- FindAllMarkers(Fresh_FRCs, assay = "SCT", slot = "counts", test.use = "MAST", latent.vars = c("S.Score", "G2M.Score", "percent.mt"), min.diff.pct = 0.2)
write.csv(Fresh_FRCs_Markers, file = "/home/zachary-wilkes/scRNAseq_Publication/Fresh_FRCs_Markers")


######

FRC.Combined.FreshVCulture <- merge(GFP_Ag_2D_LP, y = c(GFP_Ag_3D_LP, Fresh_FRCs)
                                    , add.cell.ids = c("GFPAg2DLP", "GFPAg3DLP", "FreshlyIsolated"), project = "FRC.combined.FreshVCulture")
FRC.Combined.FreshVCulture <- SetIdent(FRC.Combined.FreshVCulture, value = FRC.Combined.FreshVCulture$orig.ident)
FRC.Combined.FreshVCulture <- PrepSCTFindMarkers(FRC.Combined.FreshVCulture)
FRC.Combined.FreshVCulture <- SCTransform(FRC.Combined.FreshVCulture, verbose = T, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
FRC.Combined.FreshVCulture <- JoinLayers(FRC.Combined.FreshVCulture, assay = "RNA")

FRC.Combined.FreshVCulture <- RunPCA(FRC.Combined.FreshVCulture, assay = "SCT")

ElbowPlot(FRC.Combined.FreshVCulture)

stdv <- FRC.Combined.FreshVCulture@reductions$pca@stdev
sum.stdv <- sum(FRC.Combined.FreshVCulture@reductions$pca@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min(co1, co2)


FRC.Combined.FreshVCulture_Markers <- FindAllMarkers(FRC.Combined.FreshVCulture, assay = "SCT", slot = "data",test.use = "MAST", latent.vars = c("S.Score", "G2M.Score", "percent.mt"), min.diff.pct = 0.2)
write.csv(FRC.Combined.FreshVCulture_Markers, file = "/home/zachary-wilkes/scRNAseq_Publication/FRC.Combined.FreshVCulture_Markers.csv")

FRC.Combined.FreshVCulture <- RunUMAP(
  FRC.Combined.FreshVCulture, 
  dims = 1:15, 
  verbose = TRUE, 
  min.dist = 0.3, 
  n.epochs = 500, 
  spread     = 1, 
  n.neighbors = 100, 
  metric = "correlation", 
  assay = "SCT")


DimPlot(FRC.Combined.FreshVCulture)
# Load Necessary Libraries
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(cowplot)
library(dplyr)
library(circlize)
library(viridis)

# Define a Consistent Color Palette for Clusters Across UMAP Plots
cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                    "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# Rename 'orig.ident' for Each Dataset for Consistency
GFP_Ag_2D_LP@meta.data$orig.ident <- "2D"
GFP_Ag_3D_LP@meta.data$orig.ident <- "3D"
Fresh_FRCs@meta.data$orig.ident <- "Ex Vivo"

FRC.Combined.FreshVCulture@meta.data$orig.ident <- gsub("GFPAg2DLP", "2D", FRC.Combined.FreshVCulture@meta.data$orig.ident)
FRC.Combined.FreshVCulture@meta.data$orig.ident <- gsub("GFPAg3DLP", "3D", FRC.Combined.FreshVCulture@meta.data$orig.ident)
FRC.Combined.FreshVCulture@meta.data$orig.ident <- gsub("FreshFRCs", "Ex Vivo", FRC.Combined.FreshVCulture@meta.data$orig.ident)

# Generate UMAP Plots for Each Dataset with Enhanced Titles
umap_2D <- DimPlot(GFP_Ag_2D_LP, cols = cluster_colors) + 
  ggtitle("UMAP 2D") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

umap_3D <- DimPlot(GFP_Ag_3D_LP, cols = cluster_colors) + 
  ggtitle("UMAP 3D") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

umap_Fresh <- DimPlot(Fresh_FRCs, cols = cluster_colors) + 
  ggtitle("UMAP Ex Vivo") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

umap_Combined <- DimPlot(FRC.Combined.FreshVCulture, group.by = "orig.ident", cols = cluster_colors) + 
  ggtitle("UMAP Combined") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Define a Colorblind-Friendly Palette for Heatmaps
heatmap_colors <- viridis(3, direction = 1)

# Function to Generate Heatmaps with Legible Gene Names

generate_heatmap <- function(seurat_obj, top_n = 10, row_split_size = 10, assay_name = "SCT", 
                             test_method = "MAST", latent_vars = c("S.Score", "G2M.Score", "percent.mt"), 
                             group_by_cluster = "seurat_clusters", 
                             column_order_by = "cluster_number",
                             fresh_first = FALSE) {
  
  # 1. Identify Markers for Each Cluster
  markers <- FindAllMarkers(seurat_obj, assay = assay_name, test.use = test_method, 
                            latent.vars = latent_vars, only.pos = FALSE, 
                            min.pct = 0.40, logfc.threshold = 0.25, slot = "data")
  
  # 2. Select the Top Markers for Each Cluster
  top_markers <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = top_n, wt = avg_log2FC)
  
  # 3. Scale the Data for Selected Top Markers
  seurat_obj <- ScaleData(seurat_obj, features = top_markers$gene, assay = assay_name)
  
  # 4. Extract the Scaled Data Matrix for the Heatmap
  mat <- seurat_obj@assays[[assay_name]]$scale.data[top_markers$gene, ]
  
  # 5. Define Row Split Points for Visual Separation in the Heatmap
  num_rows <- nrow(mat)
  row_split_points <- rep(1:ceiling(num_rows / row_split_size), each = row_split_size)[1:num_rows]
  
  # 6. Reorder Columns Based on Specified Logic
  if (column_order_by == "cluster_number") {
    # Order Columns by Increasing Cluster Number
    column_order <- order(seurat_obj@meta.data[[group_by_cluster]])
  } else if (column_order_by == "fresh_2D_3D") {
    # Order Columns as Ex Vivo, 2D, 3D
    levels <- c("Ex Vivo", "2D", "3D")
    column_order <- order(factor(seurat_obj@meta.data[[group_by_cluster]], levels = levels))
  }
  
  # Apply the Column Reordering to the Data Matrix
  mat <- mat[, column_order]
  
  # 7. Order Rows Based on Grouping
  if (group_by_cluster == "orig.ident") {
    # Order Rows to Have "Ex Vivo", "2D", "3D" Clusters First
    row_order <- order(factor(top_markers$cluster, levels = c("Ex Vivo", "2D", "3D"))) 
  } else {
    # For Other Groupings, Order Based on Cluster Numbers
    cluster_order_map <- setNames(seq_along(unique(seurat_obj@meta.data[[group_by_cluster]])), 
                                  unique(seurat_obj@meta.data[[group_by_cluster]]))
    row_order <- order(cluster_order_map[top_markers$cluster])
  }
  
  # Apply the Row Reordering to the Data Matrix and Top Markers Dataframe
  mat <- mat[row_order, ]
  top_markers <- top_markers[row_order, ]
  
  # 8. Create the Heatmap Object with Enhanced Aesthetics
  heatmap_obj <- Heatmap(mat, 
                         cluster_rows = FALSE,    # Disable Row Clustering (Custom Order)
                         cluster_columns = FALSE,  # Disable Column Clustering (Custom Order)
                         show_row_names = FALSE,  # Enable Row Names (Gene Names) for Clarity # Increase Gene Name Font Size
                         show_column_names = FALSE, # Hide Column Names (Cell IDs) for Clarity
                         cluster_row_slices = TRUE, # Enable Row Clustering Within Each Split
                         cluster_column_slices = FALSE, # Disable Column Clustering Within Each Split
                         row_split = row_split_points, # Split Rows Based on Defined Points
                         column_split = seurat_obj@meta.data[[group_by_cluster]][column_order], # Split Columns by Group
                         col = heatmap_colors,     # Apply the Defined Color Palette
                         show_column_dend = FALSE,  # Hide the Column Dendrogram
                         show_row_dend = FALSE,    # Hide the Row Dendrogram
                         row_title = NULL,         # No Title for Rows
                         height = unit(5.5, "inch"),
                         
                         # Add Gene Names as Annotations on the Left
                         left_annotation = rowAnnotation(
                           genes = anno_text(rownames(mat), 
                                             location = 1,   # Center the Text
                                             just = "right",    # Align Text to the Right
                                             gp = gpar(fontsize = 10, fontface = "bold")), # Increased Font Size for Better Readability
                           annotation_width = unit(12, "mm")# Increased Width for Better Spacing
                         ),
                         
                         # Customize the Heatmap Legend
                         heatmap_legend_param = list(
                           title = "Expression", 
                           title_gp = gpar(fontsize = 12, fontface = "bold"),
                           labels_gp = gpar(fontsize = 10), 
                           legend_width = unit(8, "cm"), 
                           legend_direction = "vertical",
                           just = "center"
                         ),
                         
                         # Adjust Row Gaps for Better Visibility
                         row_gap = unit(1, "mm")      # Increase row gap for better clarity
  )
  
  # 9. Draw the Heatmap and Capture as a Grob using grid.grabExpr
  heatmap_grob <- grid.grabExpr(draw(heatmap_obj, 
                                     heatmap_legend_side = "right",    # Position the Heatmap Legend on the Right
                                     annotation_legend_side = "bottom", # Position the Annotation Legend at the Bottom
                                     padding = unit(c(10, 10, 10, 10), "mm"))) # Increase Padding for Clarity
  
  return(heatmap_grob)
}

# Generate Heatmap Grobs for Each Condition
heatmap_2D_grob <- generate_heatmap(GFP_Ag_2D_LP, top_n = 10, row_split_size = 10, 
                                    column_order_by = "cluster_number")

heatmap_3D_grob <- generate_heatmap(GFP_Ag_3D_LP, top_n = 10, row_split_size = 10,
                                    column_order_by = "cluster_number")

heatmap_Fresh_grob <- generate_heatmap(Fresh_FRCs, top_n = 5, row_split_size = 5,
                                       column_order_by = "cluster_number")

heatmap_Combined_grob <- generate_heatmap(FRC.Combined.FreshVCulture, top_n = 10, row_split_size = 10, 
                                          group_by_cluster = "orig.ident",
                                          column_order_by = "fresh_2D_3D", 
                                          fresh_first = TRUE)

# Create Title Grobs for Each Heatmap with Transparent Backgrounds
heatmap_2D_title <- ggdraw() + 
  draw_label("Heatmap of Top 10 Markers in 2D Condition", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(-19, 35, 0, 0))

heatmap_3D_title <- ggdraw() + 
  draw_label("Heatmap of Top 10 Markers in 3D Condition", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(-19, 35, 0, 0))

heatmap_Fresh_title <- ggdraw() + 
  draw_label("Heatmap of Top 5 Markers in Ex Vivo Condition", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(-19, 35, 0, 0))

heatmap_Combined_title <- ggdraw() + 
  draw_label("Heatmap of Top 10 Markers in Combined Samples", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(-19, 35, 0, 0))

# Combine Each Heatmap with Its Corresponding Title Vertically
heatmap_2D_plot <- plot_grid(heatmap_2D_title, heatmap_2D_grob, ncol = 1, rel_heights = c(0.1, 1))
heatmap_3D_plot <- plot_grid(heatmap_3D_title, heatmap_3D_grob, ncol = 1, rel_heights = c(0.1, 1))
heatmap_Fresh_plot <- plot_grid(heatmap_Fresh_title, heatmap_Fresh_grob, ncol = 1, rel_heights = c(0.1, 1))
heatmap_Combined_plot <- plot_grid(heatmap_Combined_title, heatmap_Combined_grob, ncol = 1, rel_heights = c(0.1, 1))


# Arrange All Plots in a 4x2 Grid with Consistent Spacing and Alignment
figure_1 <- plot_grid(
  plot_grid(umap_2D, heatmap_2D_plot, labels = c("A", "B"), label_size = 15, ncol = 2, align = "v", rel_widths = c(1, 1)),
  plot_grid(umap_3D, heatmap_3D_plot, labels = c("C", "D"), label_size = 15, ncol = 2, align = "v", rel_widths = c(1, 1)),
  plot_grid(umap_Fresh, heatmap_Fresh_plot, labels = c("E", "F"), label_size = 15, ncol = 2, align = "v", rel_widths = c(1, 1)),
  plot_grid(umap_Combined, heatmap_Combined_plot, labels = c("G", "H"), label_size = 15, ncol = 2, align = "v", rel_widths = c(1, 1)),
  nrow = 4, align = "v", axis = "lr",
  rel_heights = c(1, 1, 1, 1)
)

# Save the Final Figure with High Resolution Suitable for Publication
ggsave("figure_1.png", figure_1, width = 16, height = 26, dpi = 300, bg = "white")



#####
# Load necessary libraries
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(reshape2)
library(stringr)
library(ggpubr)
library(cowplot)
library(tidyr)
library(rstatix)
library(monocle3)
library(decoupleR)
library(pheatmap)
library(ggplotify)
library(viridis)
library(scales)
library(ggthemes)
library(gridExtra)
library(grid)
library(msigdbr)       # For MSigDB gene sets
library(dplyr)
library(progeny)       # For PROGENy models
library(ComplexHeatmap)
library(ggrepel)
library(EnhancedVolcano)
library(forcats)       # For reordering factors
library(ggtext)        # For improved text rendering

# -------------------------------
# Part 1: Data Preparation and Differential Expression Analysis
# -------------------------------

# Assume that FRC.Combined.FreshVCulture Seurat object is already loaded into the environment.

# Set the identity class for the Seurat object
FRC.Combined.FreshVCulture <- SetIdent(FRC.Combined.FreshVCulture, value = FRC.Combined.FreshVCulture$orig.ident)

# Define your conditions
conditions <- c("2D", "3D", "Ex Vivo")

# Prepare a list to store ranked gene lists for each comparison
ranked_gene_lists <- list()

# Perform pairwise comparisons
comparisons <- list(
  "2D_vs_ExVivo" = c("2D", "Ex Vivo"),
  "3D_vs_ExVivo" = c("3D", "Ex Vivo")
)

for (comp_name in names(comparisons)) {
  ident1 <- comparisons[[comp_name]][1]
  ident2 <- comparisons[[comp_name]][2]
  
  # Perform differential expression analysis
  markers <- FindMarkers(
    FRC.Combined.FreshVCulture,
    ident.1 = ident1,
    ident.2 = ident2,
    assay = "SCT", slot = "data",
    test.use = "MAST",
    logfc.threshold = 0,  # Include all genes
    latent.vars = c("percent.mt", "S.Score", "G2M.Score"),
    norm.method = "SCTransform",
    fc.slot = "data"
  )

  markers <- subset(markers, subset = markers$p_val_adj < 0.005)
  # Prepare ranked gene list
  markers <- markers %>%
    arrange(desc(avg_log2FC))
  
  gene_list <- markers$avg_log2FC
  names(gene_list) <- rownames(markers)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Store in the list
  ranked_gene_lists[[comp_name]] <- gene_list
}

# Check if ranked_gene_lists has data
if (length(ranked_gene_lists) == 0) {
  stop("No gene lists found for any comparisons. Cannot proceed with GSEA.")
}

# -------------------------------
# Part 2: GSEA and Visualization
# -------------------------------

# Function to filter out cell replication pathways
filter_cell_cycle_terms <- function(enrichment_result) {
  cell_cycle_terms <- c(
    
  )
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    return(enrichment_result)
  }
  filtered_result <- enrichment_result
  filtered_result@result <- enrichment_result@result[!grepl(
    paste(cell_cycle_terms, collapse = "|"),
    enrichment_result@result$Description,
    ignore.case = TRUE
  ), ]
  return(filtered_result)
}

### **GSEA GO Analysis**

# Initialize list to store GSEA results
gsea_go_results <- list()


for (comp_name in names(ranked_gene_lists)) {
  gene_list <- ranked_gene_lists[[comp_name]]
  
  # Run GSEA with gseGO
  gseGO_res <- gseGO(
    geneList = gene_list,
    ont = "BP",
    keyType = "SYMBOL",
    minGSSize = 10,
    maxGSSize = 1000,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    OrgDb = org.Mm.eg.db,
    pAdjustMethod = "BH",
    eps = 0
  )
  
  gseGO_res <- clusterProfiler::simplify(gseGO_res, cutoff = 0.4)
  
  # Filter out cell cycle-related pathways
  gsea_go_filtered <- gseGO_res
  
  # Store results
  gsea_go_results[[comp_name]] <- gsea_go_filtered
}

# Combine GSEA GO results for visualization
combine_gsea_results <- function(gsea_results_list) {
  combined_data <- do.call(rbind, lapply(names(gsea_results_list), function(comp) {
    res <- gsea_results_list[[comp]]@result
    if (!is.null(res) && nrow(res) > 0) {
      data.frame(res, Comparison = comp)
    } else {
      data.frame()
    }
  }))
  return(combined_data)
}

gsea_go_combined <- combine_gsea_results(gsea_go_results)

# Check if any results are available
if (nrow(gsea_go_combined) == 0) {
  stop("No GSEA GO results found. Cannot proceed with visualization.")
}

# Reorder the factors for better visualization
gsea_go_combined$Comparison <- factor(
  gsea_go_combined$Comparison,
  levels = c("2D_vs_ExVivo", "3D_vs_ExVivo")
)

# Create GO dotplot
go_dotplot <- ggplot(gsea_go_combined, aes(
  x = Comparison,
  y = reorder(Description, NES),
  size = -log10(p.adjust),
  color = NES
)) +
  geom_point() +
  scale_color_gradientn(colors = viridis(256)) +
  theme_bw() +
  labs(
    title = "GSEA GO Biological Process Enrichment",
    y = NULL,
    x = "Comparison",
    size = expression(-log[10](adjusted~italic(P))),
    color = "NES"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 14, face = "bold"),  # Increased x-axis label font size
    axis.title.x = element_text(size = 16, face = "bold"), # Increased x-axis title font size
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  # Adjust left margin
  theme(
    plot.margin = unit(c(1, 1, 1, 2), "cm")  # top, right, bottom, left
  ) #+
  # Wrap long y-axis labels
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
# -------------------------------
# Part 2: PCA Distance Metrics and Effect Size Calculations
# -------------------------------

# Normalize, scale, and run PCA on the Seurat object
FRC.Combined.FreshVCulture <- NormalizeData(FRC.Combined.FreshVCulture, assay = "SCT")
FRC.Combined.FreshVCulture <- ScaleData(FRC.Combined.FreshVCulture, assay = "SCT")
FRC.Combined.FreshVCulture <- RunPCA(FRC.Combined.FreshVCulture, assay = "SCT")

# Extract PCA embeddings
pca_data <- as.data.frame(FRC.Combined.FreshVCulture@reductions$pca@cell.embeddings)
pca_data <- pca_data[, 1:16]  # Use the first 16 PCs

# Add metadata
pca_data$cell_id <- rownames(pca_data)
pca_data$condition <- FRC.Combined.FreshVCulture@meta.data[rownames(pca_data), "orig.ident"]

# Calculate pairwise distances
dist_matrix <- dist(pca_data[, 1:16])
dist_df <- as.data.frame(as.table(as.matrix(dist_matrix)))
colnames(dist_df) <- c("cell1", "cell2", "distance")

# Merge condition information
dist_df <- dist_df %>%
  left_join(pca_data[, c("cell_id", "condition")], by = c("cell1" = "cell_id")) %>%
  rename(condition1 = condition) %>%
  left_join(pca_data[, c("cell_id", "condition")], by = c("cell2" = "cell_id")) %>%
  rename(condition2 = condition)

# Filter out self-comparisons and within-condition comparisons
dist_df_filtered <- dist_df %>%
  filter(cell1 != cell2, condition1 != condition2)

# Create condition pairs (sorted to avoid duplicate pairs)
condition_distances <- dist_df_filtered %>%
  rowwise() %>%
  mutate(condition_pair = paste(sort(c(condition1, condition2)), collapse = "_vs_")) %>%
  ungroup()

# Prepare data for plotting
plot_data <- condition_distances %>%
  select(condition_pair, distance)

# Ensure the condition pairs are ordered
plot_data$condition_pair <- factor(plot_data$condition_pair, levels = unique(plot_data$condition_pair))

# Calculate Effect Sizes (Cohen's d) Between Condition Pairs

# Prepare data for effect size calculation
effect_sizes <- plot_data %>%
  group_by(condition_pair) %>%
  summarize(mean_distance = mean(distance, na.rm = TRUE),
            sd_distance = sd(distance, na.rm = TRUE),
            n = n())

# Define pairwise comparisons between condition pairs
pairwise_comparisons <- combn(nrow(effect_sizes), 2, simplify = FALSE)

# Initialize list to store effect sizes
effect_sizes_list <- list()

for (comp in pairwise_comparisons) {
  group1 <- effect_sizes$condition_pair[comp[1]]
  group2 <- effect_sizes$condition_pair[comp[2]]
  
  mean1 <- effect_sizes$mean_distance[comp[1]]
  mean2 <- effect_sizes$mean_distance[comp[2]]
  
  sd1 <- effect_sizes$sd_distance[comp[1]]
  sd2 <- effect_sizes$sd_distance[comp[2]]
  
  n1 <- effect_sizes$n[comp[1]]
  n2 <- effect_sizes$n[comp[2]]
  
  # Check for sufficient data
  if(n1 < 2 | n2 < 2){
    warning(paste("Skipping comparison:", group1, "vs", group2, "- insufficient data"))
    next
  }
  
  # Calculate pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2)/(n1 + n2 - 2))
  
  # Calculate Cohen's d
  d <- (mean1 - mean2) / pooled_sd
  
  # Store the result
  effect_sizes_list[[paste(group1, group2, sep = "_vs_")]] <- data.frame(
    group1 = group1,
    group2 = group2,
    effsize = d
  )
}

# Combine effect sizes into a data frame
if(length(effect_sizes_list) > 0){
  effect_sizes_df <- bind_rows(effect_sizes_list)
  rownames(effect_sizes_df) <- NULL  # Clean row names
} else {
  effect_sizes_df <- data.frame()
  warning("No effect sizes were calculated due to insufficient data in all comparisons.")
}

# Visualization with Violin Plots and Mean Distance Annotations


# Create the violin plot
pca_distance_plot <- ggplot(plot_data, aes(
  x = condition_pair,
  y = distance,
  fill = NULL
)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  # Add mean distance annotations
  geom_text(
    data = effect_sizes,
    aes(
      x = condition_pair,
      y = max(plot_data$distance, na.rm = TRUE) * 1.05,
      label = paste0("Mean = ", round(mean_distance, 2))
    ),
    vjust = 0,
    size = 5,  # Increased annotation text size
    inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(
    title = "PCA Distances Between Conditions",
    x = "Condition Pair",
    y = "Distance"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) 

# Adjust the PCA distance plot to be wider using ggdraw
library(cowplot)

pca_distance_plot_wide <- ggdraw() +
  draw_plot(pca_distance_plot, x = -0.6, y = -0.1, width = 1.65, height = 1.15)
# Display effect sizes between specific condition pairs
target_comparisons <- c("2D_vs_3D", "2D_vs_Ex Vivo", "3D_vs_Ex Vivo")

effect_sizes_to_display <- effect_sizes_df %>%
  filter(paste(group1, group2, sep = "_vs_") %in% target_comparisons |
           paste(group2, group1, sep = "_vs_") %in% target_comparisons)

# Adjust labels for both group1_vs_group2 and group2_vs_group1
effect_sizes_to_display <- effect_sizes_to_display %>%
  mutate(condition_pair_label = ifelse(
    paste(group1, group2, sep = "_vs_") %in% target_comparisons,
    paste(group1, group2, sep = "_vs_"),
    paste(group2, group1, sep = "_vs_")
  ))

# Add effect size annotations if any are available
if(nrow(effect_sizes_to_display) > 0){
  # Determine y positions for annotations
  annotation_y_position <- max(plot_data$distance, na.rm = TRUE) * 1.1
  
  # Create a mapping from condition_pair to effect size label
  effect_size_mapping <- effect_sizes_to_display %>%
    mutate(
      label = paste0("d = ", round(effsize, 2))
    )
  
  # Add effect size annotations
  pca_distance_plot_wide <- pca_distance_plot_wide +
    geom_text(
      data = effect_size_mapping,
      aes(
        x = condition_pair_label,
        y = annotation_y_position,
        label = label
      ),
      vjust = 0,
      size = 4,
      inherit.aes = FALSE
    )
  
  # Adjust y-axis limits to accommodate the annotations
  pca_distance_plot_wide <- pca_distance_plot_wide +
    ylim(NA, max(plot_data$distance, na.rm = TRUE) * 1.2)
}


  
#pca_distance_plot_wide <- pca_distance_plot_wide + theme(
    #plot.margin = unit(c(5, 5, 1, 1), "cm")) #+
# -------------------------------
# Part 3: Trajectory Analysis with Monocle 3
# -------------------------------

# Convert Seurat object to Monocle 3 object
FRC_monocle <- as.cell_data_set(FRC.Combined.FreshVCulture)

# Transfer cluster information
FRC_monocle@clusters$UMAP$clusters <- FRC.Combined.FreshVCulture@meta.data$seurat_clusters
FRC_monocle@colData$seurat_clusters <- FRC.Combined.FreshVCulture@meta.data$seurat_clusters
FRC_monocle@colData$orig.ident <- FRC.Combined.FreshVCulture@meta.data$orig.ident

# Learn the trajectory graph
FRC_monocle <- cluster_cells(FRC_monocle, reduction_method = "UMAP")
FRC_monocle <- learn_graph(FRC_monocle, use_partition = TRUE)

# Order cells with root cells (e.g., "Ex Vivo")
root_cells <- colnames(FRC_monocle)[FRC_monocle@colData$orig.ident == "Ex Vivo"]

# Verify root cells
if(length(root_cells) == 0){
  stop("No root cells found with orig.ident == 'Ex Vivo'. Please check the identifier.")
}

FRC_monocle <- order_cells(FRC_monocle, root_cells = root_cells)

trajectory_plot_pseudotime <- plot_cells(
FRC_monocle,
color_cells_by = "pseudotime",
label_groups_by_cluster = FALSE,
label_leaves = T,
label_branch_points = T,
cell_size = 0.5,
graph_label_size = 3
) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  ggtitle("Trajectory Colored by Pseudotime") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),           # Increased axis text size
    axis.title = element_text(size = 16, face = "bold"),          # Increased axis title size
    legend.title = element_text(size = 14, face = "bold"),        # Increased legend title size
    legend.text = element_text(size = 12)                         # Increased legend text size
  ) +
  labs(color = "Pseudotime")  # Changed legend title

# Plot the trajectory colored by condition (orig.ident)
# Load necessary libraries
library(ggplot2)
library(cowplot)
library(monocle3)
library(ggrepel)
library(dplyr)

# Create the trajectory plot without group labels
trajectory_plot_condition <- plot_cells(
  FRC_monocle,
  color_cells_by = "orig.ident",
  label_cell_groups = FALSE,   # Do not label cell groups here
  label_groups_by_cluster = FALSE,
  label_leaves = T,
  label_branch_points = T,
  cell_size = 0.5,
  graph_label_size = 3
) +
  theme_minimal() +
  ggtitle("Trajectory Colored by Sample ID") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.size = unit(1.5, 'lines')
  ) +
  guides(
    color = guide_legend(
      title = "Sample ID",
      override.aes = list(size = 5)
    )
  )

# Extract the embedding coordinates and metadata
embedding_df <- as.data.frame(reducedDims(FRC_monocle)$UMAP)
colnames(embedding_df) <- c("UMAP_1", "UMAP_2")
embedding_df$orig.ident <- colData(FRC_monocle)$orig.ident

# Calculate the centers of each group (Sample ID)
label_df <- embedding_df %>%
  group_by(orig.ident) %>%
  summarize(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

# Add the labels to the plot with bold font
trajectory_plot_condition <- trajectory_plot_condition +
  geom_text_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = orig.ident),
    size = 6,
    fontface = "bold",
    segment.color = NA    # Remove line segments
  )

trajectory_plot_condition <- trajectory_plot_condition +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),           # Increased axis text size
    axis.title = element_text(size = 16, face = "bold"),          # Increased axis title size
    legend.title = element_text(size = 14, face = "bold"),        # Increased legend title size
    legend.text = element_text(size = 12)                         # Increased legend text size
  )
# -------------------------------
# Part 4: MLM and VIPER Analysis
# -------------------------------

# Ensure that the Seurat object is normalized
FRC.Combined.FreshVCulture <- NormalizeData(FRC.Combined.FreshVCulture, assay = "RNA")

# Get the PROGENy model for mouse
net <- progeny::model_mouse_full
net$gene <- as.character(net$gene)

# Extract the expression matrix
mat <- as.matrix(FRC.Combined.FreshVCulture@assays$SCT$data)

# Run MLM analysis
acts <- run_mlm(
  mat = mat,
  network = net,
  .source = 'pathway',
  .target = 'gene',
  .mor = 'weight',
  minsize = 0,
  na.rm = TRUE
)

# Add pathway activities as a new assay in Seurat object
FRC.Combined.FreshVCulture[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  t() %>%
  Seurat::CreateAssayObject()

# Change default assay to pathwaysmlm
DefaultAssay(FRC.Combined.FreshVCulture) <- "pathwaysmlm"

# Scale the data
FRC.Combined.FreshVCulture <- ScaleData(FRC.Combined.FreshVCulture)

# Extract activities and prepare data for heatmap
df_mlm <- GetAssayData(FRC.Combined.FreshVCulture, slot = "scale.data") %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(condition = FRC.Combined.FreshVCulture@meta.data$orig.ident) %>%
  group_by(condition) %>%
  summarise(across(everything(), mean))

# Transform to matrix
mlm_mat <- df_mlm %>%
  column_to_rownames("condition") %>%
  as.matrix()

# Load ComplexHeatmap for advanced heatmap plotting
library(ComplexHeatmap)

mlm_heatmap <- Heatmap(
  mlm_mat,
  name = "Expression",  # Set legend title to "Expression"
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = viridis(50),
  border = NA,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),  # Increased legend title font size
    labels_gp = gpar(fontsize = 10),                     # Increased legend labels font size             # Ensure title is centered above the legend
    title_gap = unit(20, "mm")
  )
)

# Convert MLM heatmap to grob
mlm_heatmap_grob <- grid.grabExpr(draw(mlm_heatmap))

# Create title grob for MLM heatmap
mlm_heatmap_title <- ggdraw() +
  draw_label("MLM Pathway Activity Heatmap", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_void()

# Combine MLM heatmap with title
mlm_heatmap_plot <- plot_grid(
  mlm_heatmap_title,
  mlm_heatmap_grob,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# VIPER Analysis
net_tf <- decoupleR::get_collectri(organism = 'mouse')

# Run VIPER analysis
viper_results <- run_viper(
  mat = mat,
  network = net_tf,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 0,
  pleiotropy = TRUE,
  cores = 4  # Adjust cores based on your system
)

# Add VIPER activities as a new assay in Seurat object
FRC.Combined.FreshVCulture[['viper']] <- viper_results %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  t() %>%
  Seurat::CreateAssayObject()

# Change default assay to viper
DefaultAssay(FRC.Combined.FreshVCulture) <- "viper"

# Scale the data
FRC.Combined.FreshVCulture <- ScaleData(FRC.Combined.FreshVCulture)

# Extract activities and prepare data for heatmap
df_viper <- GetAssayData(FRC.Combined.FreshVCulture, slot = "scale.data") %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(condition = FRC.Combined.FreshVCulture@meta.data$orig.ident) %>%
  group_by(condition) %>%
  summarise(across(everything(), mean))

# Get top variable transcription factors
top_n_tfs <- 50
tfs_sd <- apply(df_viper[,-1], 2, sd)
top_tfs <- names(sort(tfs_sd, decreasing = TRUE))[1:top_n_tfs]

# Subset data to top transcription factors
viper_mat <- df_viper %>%
  select(condition, all_of(top_tfs)) %>%
  column_to_rownames("condition") %>%
  as.matrix()

viper_heatmap <- Heatmap(
  viper_mat,
  name = "Expression",  # Set legend title to "Expression"
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = viridis(50),
  border = NA,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),  # Increased legend title font size
    labels_gp = gpar(fontsize = 10),# Increased legend labels font size             # Ensure title is centered above the legend
    title_gap = unit(20, "mm")
  )
)

# Convert VIPER heatmap to grob
viper_heatmap_grob <- grid.grabExpr(draw(viper_heatmap))

# Create title grob for VIPER heatmap
viper_heatmap_title <- ggdraw() +
  draw_label("VIPER Transcription Factor Activity Heatmap", fontface = 'bold', size = 14, hjust = 0.5) +
  theme_void()

# Combine VIPER heatmap with title
viper_heatmap_plot <- plot_grid(
  viper_heatmap_title,
  viper_heatmap_grob,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# -------------------------------
# Part 5: Combine Plots into Figure 2
# -------------------------------
go_dotplot <- go_dotplot +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

pca_distance_plot <- pca_distance_plot_wide + theme(legend.position = "none")

trajectory_plot_pseudotime <- trajectory_plot_pseudotime + theme(legend.position = "right")

trajectory_plot_condition <- trajectory_plot_condition + theme(legend.position = "right")

# First row: GO dotplot and PCA distance plot
dotplot_row <- plot_grid(
  go_dotplot,
  pca_distance_plot,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1),  # Allocate more width to the dotplot
  label_size = 24
)

# Second row: MLM and VIPER plots in two columns
mlm_viper_row <- plot_grid(
  mlm_heatmap_plot,
  viper_heatmap_plot,
  labels = c("C", "D"),
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1),
  label_size = 24
)

# Third row: Trajectory plots in two columns
trajectory_row <- plot_grid(
  trajectory_plot_pseudotime,
  trajectory_plot_condition,
  labels = c("E", "F"),
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1),
  label_size = 24
)

# Combine all rows into the final figure
figure_2_combined <- plot_grid(
  dotplot_row,
  mlm_viper_row,
  trajectory_row,
  ncol = 1,
  rel_heights = c(1, 1, 1)
)

# Save the combined figure with increased width
ggsave(
  "Figure_2_Combined.png",
  figure_2_combined,
  width = 25,  # Increased width
  height = 25,
  dpi = 300,
  bg = "white"
)
# -------------------------------
# Additional Code: Volcano Plots
# -------------------------------

# Load necessary libraries
library(Seurat)
library(EnhancedVolcano)
library(ggplot2)
library(cowplot)

# -------------------------------
# Part 1: Differential Expression Analysis for Volcano Plots
# -------------------------------

# Prepare a list to store markers for volcano plots
volcano_markers <- list()

# Perform pairwise comparisons
comparisons <- list(
  "A: 2D vs Ex Vivo" = c("2D", "Ex Vivo"),
  "B: 3D vs Ex Vivo" = c("3D", "Ex Vivo"),
  "C: 2D vs 3D" = c("2D", "3D")
)

for (comp_name in names(comparisons)) {
  ident1 <- comparisons[[comp_name]][1]
  ident2 <- comparisons[[comp_name]][2]
  
  # Perform differential expression analysis
  markers <- FindMarkers(
    FRC.Combined.FreshVCulture,
    ident.1 = ident1,
    ident.2 = ident2,
    slot = "counts",
    assay = "SCT",
    test.use = "MAST",
    min.pct = 0.1,            # Set a minimum percentage
    logfc.threshold = 0.25,   # Set a logFC threshold
    latent.vars = c("percent.mt", "S.Score", "G2M.Score")
  )
  
  volcano_markers[[comp_name]] <- markers
}

# -------------------------------
# Part 2: Generate and Customize Volcano Plots
# -------------------------------

# Function to create customized volcano plots
create_volcano_plot <- function(data, title) {
  # Identify top genes based on adjusted p-value and log2 fold change
  
  data$col <- "grey"
  
  # Define significance criteria
  pcut <- 1e-05
  fccut <- 1.0
  
  pos_genes <- which(data$p_val_adj < pcut & data$avg_log2FC >= fccut)
  neg_genes <- which(data$p_val_adj < pcut & data$avg_log2FC <= -fccut)
  pval_only_genes <- which(data$padj < pcut & abs(data$avg_log2FC) < fccut)
  
  data$col[pval_only_genes] <- "blue"
  # Negative significant: red
  data$col[neg_genes] <- "red"
  # Positive significant: green
  data$col[pos_genes] <- "green"
  
  custom_col <- data$col
  
  names(custom_col)[custom_col == 'grey'] <- 'NS'
  names(custom_col)[custom_col == 'blue'] <- 'P-Value Significant'
  names(custom_col)[custom_col == 'red'] <- 'Significant and Negative'
  names(custom_col)[custom_col == 'green'] <- 'Significant and Positive'
  
  top_up <- data %>%
    mutate(gene = rownames(data)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(20) %>%
    pull(gene)
  
  top_down <- data %>%
    mutate(gene = rownames(data)) %>%
    arrange(avg_log2FC) %>%
    head(20) %>%
    pull(gene)
  
  top_genes <- c(top_up, top_down)
  
  
  
  EnhancedVolcano(
    data,
    lab = rownames(data),
    selectLab = top_genes,  # Label only the selected top genes
    x = "avg_log2FC",
    y = "p_val_adj",
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~ adjusted~italic(P)),
    title = title,
    arrowheads = FALSE,
    pCutoff = pcut,
    FCcutoff = fccut,
    pointSize = 1,
    subtitle = NULL,
    # Remove xlim and ylim first to let EnhancedVolcano autoscale
    xlim = c(-15, 15),
    ylim = c(0, 326),
    labSize = 3,
    colCustom = custom_col,
    colAlpha = 0.6,
    drawConnectors = TRUE,
    max.overlaps = Inf,# Enable connectors to better show labels
    widthConnectors = 0.5,
    colConnectors = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = "full",
    borderWidth = 0.5,
    borderColour = "black",
    caption = NULL
  ) +
    theme_bw(base_size = 14) +  # Use a clean theme
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )
}

# Generate volcano plots for each comparison
volcano_plots <- list()
for (comp_name in names(volcano_markers)) {
  plot <- create_volcano_plot(
    data = volcano_markers[[comp_name]],
    title = comp_name
  )
  volcano_plots[[comp_name]] <- plot
}

# -------------------------------
# Part 3: Combine Volcano Plots into Figure 3
# -------------------------------

# Arrange the volcano plots in a grid

figure_3 <- plot_grid(
  volcano_plots[["A: 2D vs Ex Vivo"]],
  volcano_plots[["B: 3D vs Ex Vivo"]],
  volcano_plots[["C: 2D vs 3D"]],
  labels = c("A", "B", "C"),
  ncol = 1,
  align = "hv",
  axis = "tblr",
  rel_heights = c(1, 1, 1)
)

# Save the combined figure with increased width
ggsave(
  filename = "Figure_3_Volcano_Plots.png",
  plot = figure_3,
  width = 10,   # Increased width
  height = 18,
  dpi = 300,
  bg = "white"
)
