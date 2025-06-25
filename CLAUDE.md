# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a comprehensive bioinformatics pipeline focused on RNA-seq analysis (both bulk and single-cell), with support for proteomics and general genomics workflows. The pipeline supports multiple workflow managers (Nextflow, Snakemake) and execution environments (local, HPC, cloud).

### Current Focus: LNSC Bayesian Hyperparameter Optimization

The project currently includes GPU-accelerated Bayesian hyperparameter optimization for single-cell RNA sequencing (scRNA-seq) clustering of mouse lymph node stromal cells (LNSCs). This work aims to find optimal clustering parameters using validated biological markers as ground truth.

**Key Details:**
- **Data**: Mouse lymph node stromal cells (Pdpn+ CD31- CD45- sorted populations)
- **Samples**: Three conditions - 2D in vitro, 3D in vitro, and in vivo
- **GPU**: NVIDIA RTX 4090 with CUDA
- **Goal**: Optimize clustering to identify subtle cell populations with gradient expression patterns

## Key Commands

### Environment Setup
```bash
# Full environment setup (conda, containers, test data)
./bin/setup.sh

# Setup only conda environments
./bin/setup.sh --conda-only

# Check dependencies
./bin/setup.sh --check-deps

# Setup specific environments
./bin/setup.sh --envs-only
```

### Running Analysis
```bash
# Run default Nextflow pipeline
./bin/run_analysis.sh

# Run with specific workflow and profile
./bin/run_analysis.sh -w nextflow -p docker -i data/samples.tsv

# Run Snakemake workflow
./bin/run_analysis.sh -w snakemake -p conda

# Dry run (test without execution)
./bin/run_analysis.sh --dry-run

# Run with custom resources
./bin/run_analysis.sh --cores 8 --memory 16GB
```

### Direct Workflow Execution
```bash
# Nextflow
nextflow run workflows/nextflow/main.nf --input data/samples.tsv --genome data/reference/genome.fasta -profile docker

# Snakemake
snakemake --snakefile workflows/snakemake/Snakefile --cores 4 --use-conda
```

### Cleanup
```bash
# Clean temporary files
./bin/clean.sh

# Clean logs and temp files
./bin/clean.sh -l -t

# Show what would be cleaned
./bin/clean.sh --dry-run

# Clean everything (careful!)
./bin/clean.sh --all
```

## Architecture and Structure

### Directory Organization
- `bin/` - Operational scripts (setup, run, clean)
- `code/` - Source code organized by language:
  - `R/` - R analysis scripts and modules
  - `python/` - Python utilities
  - `shell/` - Shell utilities
- `config/` - YAML configuration files
- `containers/` - Docker/Singularity definitions
- `data/` - Input data organization (raw, processed, references)
- `envs/` - Environment specifications (conda, renv)
- `workflows/` - Nextflow and Snakemake workflows
- `results/` - Analysis outputs
- `logs/` - Execution logs

### Key Analysis Modules

**Single-cell RNA-seq** (`code/R/main.R`):
- Orchestrates Seurat-based analysis
- Modular structure with source files for each step
- Default QC: min_genes=200, max_genes=6000, max_mt=20%

**Bulk RNA-seq** (`code/R/src/Bulk_RNA_Seq_Git.R`):
- DESeq2/edgeR analysis
- Comprehensive visualization
- Publication-ready outputs

**Workflow Components**:
1. Quality Control (FastQC)
2. Read Trimming (Trimmomatic)
3. Alignment (STAR, HISAT2)
4. Quantification (Salmon, Kallisto)
5. Differential Expression (DESeq2, edgeR)
6. Report Generation (MultiQC)

### Configuration System

**Main Configuration** (`config/project.yml`):
- Project metadata and versioning
- Default workflow settings
- Tool selections
- Environment configurations

**Analysis Parameters** (`config/params.yml`):
- QC thresholds
- Tool-specific parameters
- Analysis settings

**Resource Configuration** (`config/resources.yml`):
- CPU, memory, and time allocations
- Environment-specific settings (local, HPC, cloud)

### Execution Profiles

Available profiles for different environments:
- `docker` - Docker containers
- `singularity` - Singularity containers
- `podman` - Podman containers
- `conda` - Conda/Mamba environments
- `hpc` - HPC cluster execution
- `cloud` - Cloud platform execution
- `test` - Test profile with minimal resources
- `debug` - Debugging with verbose output

### Important Development Notes

1. **Modular R Code**: The R analysis code is modularized with separate source files for different analysis steps (clustering, differential expression, visualization, etc.)

2. **Configuration-Driven**: Most analysis parameters are defined in YAML files rather than hardcoded

3. **Multi-Environment Support**: The pipeline can run with different package managers and containers

4. **Workflow Flexibility**: Both Nextflow and Snakemake are fully supported with equivalent functionality

5. **Resource Management**: Dynamic resource allocation based on process requirements and available infrastructure

6. **Sample Organization**: Input samples should be specified in TSV format with paths to FASTQ files

7. **Reference Data**: Reference genomes, indices, and databases are configured in `config/project.yml`

## LNSC Optimization Specific Guidelines

### Technical Stack for GPU Optimization
- **Python 3.10+** with RAPIDS ecosystem
- **Key Libraries**: `rapids-singlecell`, `cupy`, `cudf`, `cuml`, `scanpy`, `optuna`
- **GPU Memory**: Pre-allocate with RMM for 16GB pool

### Biological Markers (Ground Truth)
```python
# Tier 1 - Binary Markers (High Specificity)
tier1_markers = {
    'FDC': ['Cr2', 'Coch', 'Fcgr2b'],
    'MRC': ['Madcam1', 'Tnfsf11'],  # RANKL
    'Medullary': ['Inmt'],
    'CD34_stromal': ['Cd34', 'Pi16', 'Cxcl14'],
    'TRC_core': ['Ccl19', 'Ccl21', 'Il7', 'Pdpn']
}

# Tier 2 - Gradient Markers
gradient_markers = {
    'TRC_subtypes': {'Ccl19': {'high': >2.5, 'low': 1.0-2.5}},
    'Inflammatory': ['Cxcl9', 'Cxcl10'],
    'Activation': ['Nr4a1']
}
```

### Critical Data Handling Notes
1. **Data arrives pre-processed from Seurat** - Do NOT re-normalize
2. **H5AD format contains**: Normalized/scaled expression, PCA embeddings, cell metadata
3. **Gene names**: May need capitalization matching between R/Python
4. **Expression levels**: Are log-normalized from Seurat

### GPU Optimization Workflow
```python
# Memory management
import rmm
rmm.reinitialize(pool_allocator=True, initial_pool_size=2**30 * 16)  # 16GB

# Data loading (already normalized from Seurat)
adata = sc.read_h5ad(path)  # Do NOT normalize again
X_gpu = cp.array(adata.X.toarray() if sparse else adata.X)

# Clear GPU memory between trials
cp.get_default_memory_pool().free_all_blocks()
```

### Success Criteria for Optimization
- Clusters show high marker specificity (>2 fold enrichment)
- CCL19hi/lo populations are distinguished
- Rare populations (Nr4a1+) are detected
- Parameters transfer well to in vitro samples

### Common Issues and Solutions
- **OOM errors**: Use batch processing or reduce PCA dimensions
- **Marker mismatches**: Check gene presence with `if marker in adata.var_names:`
- **GPU verification**: `cp.cuda.runtime.getDeviceCount()`
- **10x dropout**: Absence doesn't mean non-expression