"""
LNSC Bayesian Hyperparameter Optimization Module
GPU-accelerated clustering optimization for single-cell RNA-seq data
"""

from .bayesian_optimizer import LNSCOptimizer
from .gpu_clustering import GPUClusterer
from .marker_scoring import MarkerScorer
from .utils import (
    setup_gpu_memory,
    clear_gpu_memory,
    load_preprocessed_data,
    save_results,
    load_config,
    get_gpu_memory_info,
    monitor_gpu_memory,
    validate_markers
)

__version__ = "1.0.0"
__all__ = [
    "LNSCOptimizer",
    "GPUClusterer",
    "MarkerScorer",
    "setup_gpu_memory",
    "clear_gpu_memory",
    "load_preprocessed_data",
    "save_results",
    "load_config",
    "get_gpu_memory_info",
    "monitor_gpu_memory",
    "validate_markers"
]