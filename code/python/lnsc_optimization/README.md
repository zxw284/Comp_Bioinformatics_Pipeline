# LNSC Bayesian Hyperparameter Optimization

GPU-accelerated Bayesian optimization for single-cell RNA-seq clustering of mouse lymph node stromal cells (LNSCs).

## Overview

This module implements a biologically-driven optimization approach to find optimal clustering parameters for identifying subtle cell populations in LNSC data. It uses validated markers from Rodda et al. 2018 as ground truth to guide the optimization process.

## Features

- **GPU Acceleration**: Uses RAPIDS ecosystem for fast clustering on NVIDIA GPUs
- **Biological Validation**: Marker-based scoring using established LNSC markers
- **Bayesian Optimization**: Efficient parameter search using Optuna
- **Multi-sample Support**: Validate optimized parameters across conditions
- **Memory Efficient**: Smart GPU memory management for large datasets

## Requirements

- NVIDIA GPU with CUDA support (tested on RTX 4090)
- CUDA 11.5 or higher
- Python 3.10+
- See `envs/conda/lnsc_gpu.yml` for full dependencies

## Installation

1. Create conda environment:
```bash
conda env create -f envs/conda/lnsc_gpu.yml
conda activate lnsc_gpu
```

2. Verify GPU setup:
```bash
python -c "import cupy; print(f'GPUs: {cupy.cuda.runtime.getDeviceCount()}')"
```

## Usage

### Quick Start

```bash
# Run optimization on preprocessed data
./bin/run_lnsc_optimization.sh data/processed/lnsc_invivo.h5ad

# With validation samples
./bin/run_lnsc_optimization.sh \
    -v "data/processed/lnsc_2d.h5ad data/processed/lnsc_3d.h5ad" \
    data/processed/lnsc_invivo.h5ad
```

### Python API

```python
from lnsc_optimization import LNSCOptimizer

# Initialize optimizer
optimizer = LNSCOptimizer('config/lnsc_optimization.json')

# Run optimization
results = optimizer.optimize(n_trials=100)

# Validate on new samples
validation = optimizer.validate_on_samples(['sample1.h5ad', 'sample2.h5ad'])
```

### Configuration

Edit `config/lnsc_optimization.json` to customize:
- Marker sets (tier1, gradient)
- Parameter search ranges
- Scoring weights
- GPU memory allocation

## Input Data Format

Data must be preprocessed h5ad files from Seurat containing:
- Normalized and scaled expression matrix (`adata.X`)
- PCA embeddings (`adata.obsm['X_pca']`)
- Cell metadata (`adata.obs`)

**Important**: Do NOT re-normalize the data - it arrives already processed from Seurat.

## Outputs

1. **Optimization Results** (`*_results.json`):
   - Best parameters (n_pcs, n_neighbors, min_dist, resolution)
   - Optimization scores
   - Trial history

2. **Clustered Data** (`*_clustered.h5ad`):
   - Original data with optimized clustering
   - UMAP coordinates
   - Leiden cluster assignments

3. **Reports**:
   - Summary report (`*_summary.txt`)
   - Trial details (`*_trials.csv`)
   - Marker expression reports

## Biological Markers

### Tier 1 (Binary, High Specificity)
- **FDC**: Cr2, Coch, Fcgr2b
- **MRC**: Madcam1, Tnfsf11 (RANKL)
- **Medullary**: Inmt
- **CD34+ stromal**: Cd34, Pi16, Cxcl14
- **TRC core**: Ccl19, Ccl21, Il7, Pdpn

### Tier 2 (Gradient Expression)
- **TRC subtypes**: Ccl19 (high >2.5, low 1.0-2.5)
- **Inflammatory**: Cxcl9, Cxcl10
- **Activation**: Nr4a1

## Troubleshooting

### GPU Memory Issues
```python
# Reduce batch size or PCA dimensions
config['optimization']['param_ranges']['n_pcs']['max'] = 30

# Clear GPU memory between trials
import cupy as cp
cp.get_default_memory_pool().free_all_blocks()
```

### Missing Markers
- Check gene name capitalization (mouse vs human)
- Some genes may have synonyms (e.g., Tnfsf11 = RANKL)
- 10x platform has dropout - absence doesn't mean non-expression

### Performance Tips
- Start with 50-100 trials for exploration
- Use n_jobs=1 for GPU (parallel trials need more memory)
- Monitor GPU usage with `nvidia-smi -l 1`

## References

- Rodda et al. (2018) - Single-cell RNA sequencing of lymph node stromal cells
- RAPIDS.ai - GPU-accelerated data science
- Optuna - Hyperparameter optimization framework