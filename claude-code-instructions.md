# Claude Code Instructions: LNSC Bayesian Hyperparameter Optimization Project

## Project Overview
You are assisting with a GPU-accelerated Bayesian hyperparameter optimization project for single-cell RNA sequencing (scRNA-seq) clustering of mouse lymph node stromal cells (LNSCs). The goal is to find optimal clustering parameters using validated biological markers as ground truth.

## Biological Context
- **Data**: Mouse lymph node stromal cells (Pdpn+ CD31- CD45- sorted populations)
- **Platform**: 10x Genomics Chromium (same as reference study)
- **Samples**: Three conditions - 2D in vitro, 3D in vitro, and in vivo
- **Challenge**: Identifying subtle cell populations with gradient expression patterns
- **Preprocessing Status**: Data arrives as processed Seurat objects with:
  - QC completed (mitochondrial/complexity filtering, doublet removal)
  - Normalized and scaled (likely SCTransform)
  - Initial PCA computed (50 dims)
  - Ready for downstream analysis

## Technical Stack
- **GPU**: NVIDIA RTX 4090 with CUDA
- **Language**: Python 3.10+
- **Key Libraries**:
  - `rapids-singlecell` (GPU acceleration)
  - `cupy`, `cudf`, `cuml` (RAPIDS ecosystem)
  - `scanpy` (single-cell analysis)
  - `optuna` (Bayesian optimization)
  - `torch` (if needed for deep learning approaches)

## Ground Truth Markers
Based on validated Rodda et al. 2018 markers (extensively replicated 2018-2025):

### Tier 1 - Binary Markers (High Specificity)
```python
tier1_markers = {
    'FDC': ['Cr2', 'Coch', 'Fcgr2b'],
    'MRC': ['Madcam1', 'Tnfsf11'],  # RANKL
    'Medullary': ['Inmt'],
    'CD34_stromal': ['Cd34', 'Pi16', 'Cxcl14'],
    'TRC_core': ['Ccl19', 'Ccl21', 'Il7', 'Pdpn']
}
```

### Tier 2 - Gradient Markers
```python
gradient_markers = {
    'TRC_subtypes': {'Ccl19': {'high': >2.5, 'low': 1.0-2.5}},
    'Inflammatory': ['Cxcl9', 'Cxcl10'],
    'Activation': ['Nr4a1']
}
```

## Primary Objectives
1. Optimize clustering parameters (PCA dims, UMAP settings, Leiden resolution)
2. Maximize marker-based cluster validation scores
3. Compare core (20 genes) vs extended (50+ genes) marker sets
4. Transfer optimized parameters to 2D/3D in vitro samples

## Code Architecture

### 1. Data Loading and GPU Preparation
```python
# Data arrives from R/Seurat pipeline as h5ad
# Seurat -> h5ad conversion preserves:
# - Normalized/scaled expression matrix
# - PCA embeddings (reductions)
# - Cell metadata
# - Feature information

adata_cpu = sc.read_h5ad(path)  # Already normalized/scaled
X_gpu = cp.array(adata.X.toarray() if sparse else adata.X)

# Note: adata.X contains scaled data
# Raw counts may be in adata.raw if needed
# PCA may already be in adata.obsm['X_pca']
```

### 2. Optimization Pipeline
- Use RAPIDS for GPU acceleration
- Implement marker-based scoring functions
- Run Optuna trials with pruning
- Track both biological and technical metrics

### 3. Evaluation Metrics
- Marker specificity scores
- Gradient separation quality
- Cluster stability
- (Optional) ARI if annotations available

## Key Technical Considerations

### Memory Management
```python
# Pre-allocate GPU memory
rmm.reinitialize(pool_allocator=True, initial_pool_size=2**30 * 16)  # 16GB

# Clear between trials
cp.get_default_memory_pool().free_all_blocks()
```

### Batch Processing
- Process large datasets in chunks if needed
- Use `with cp.cuda.Device(0):` for explicit GPU control

### Error Handling
- Common issues: OOM errors, sparse matrix conversions, gene name mismatches
- Always check marker presence: `if marker in adata.var_names:`
- Implement fallbacks for missing markers

## Troubleshooting Guide

### GPU Issues
- Verify CUDA with: `cp.cuda.runtime.getDeviceCount()`
- Check memory: `cp.get_default_memory_pool().used_bytes()`
- Fallback to CPU for specific operations if needed

### Biological Validation
- Some markers may have synonyms (e.g., Tnfsf11 = RANKL)
- Expression levels are log-normalized (from Seurat)
- 10x platform has dropout - absence doesn't mean non-expression
- Gene names may need capitalization matching between R/Python

### Optimization Tips
- Start with 50-100 trials for parameter exploration
- Use pruning to stop bad trials early
- Save intermediate results: `study.trials_dataframe()`

## Example Debugging Session
```python
# Check marker expression
sc.pl.dotplot(adata, tier1_markers['FDC'], groupby='leiden')

# Verify GPU acceleration
import time
start = time.time()
# operation
print(f"GPU time: {time.time() - start:.2f}s")

# Monitor optimization
from optuna.visualization import plot_optimization_history
plot_optimization_history(study)
```

## Success Criteria
- Identified clusters show high marker specificity (>2 fold enrichment)
- CCL19hi/lo populations are distinguished
- Rare populations (Nr4a1+) are detected
- Parameters transfer well to in vitro samples

## Additional Context
- User has strong R/Seurat background
- No wet lab validation planned (relying on published validations)
- Goal is to solve challenging clustering problems in novel in vitro systems
- **CRITICAL**: Data arrives pre-processed from Seurat (QC, normalized, scaled) - do NOT re-normalize

## R-Python Workflow Notes
```r
# Seurat export (already done by user):
SaveH5Seurat(seurat_obj, "data.h5Seurat")
Convert("data.h5Seurat", dest = "h5ad")
```

```python
# Python receives:
adata = sc.read_h5ad("data.h5ad")
# adata.X is already normalized/scaled
# adata.obsm['X_pca'] has Seurat's PCA
```

Remember: The biological insight is as important as the technical implementation. Always validate that the optimized parameters make biological sense!