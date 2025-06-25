#!/usr/bin/env python3
"""
Utility functions for GPU memory management and data handling
"""

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Union

import cupy as cp
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

try:
    import rmm
    HAS_RMM = True
except ImportError:
    HAS_RMM = False
    logging.warning("RMM not available. GPU memory pool management will be limited.")

logger = logging.getLogger(__name__)


def setup_gpu_memory(pool_size_gb: float = 16.0, device_id: int = 0):
    """
    Setup GPU memory pool for efficient allocation
    
    Parameters:
    -----------
    pool_size_gb : float
        Size of memory pool in GB
    device_id : int
        GPU device ID to use
    """
    # Set device
    cp.cuda.Device(device_id).use()
    
    if HAS_RMM:
        # Initialize RMM memory pool
        pool_size_bytes = int(pool_size_gb * 1024**3)
        rmm.reinitialize(
            pool_allocator=True,
            initial_pool_size=pool_size_bytes,
            devices=[device_id]
        )
        logger.info(f"Initialized RMM memory pool with {pool_size_gb} GB on device {device_id}")
    else:
        # Use CuPy's memory pool
        mempool = cp.get_default_memory_pool()
        mempool.set_limit(size=int(pool_size_gb * 1024**3))
        logger.info(f"Set CuPy memory pool limit to {pool_size_gb} GB on device {device_id}")
    
    # Log GPU info
    device = cp.cuda.Device(device_id)
    free_mem, total_mem = device.mem_info
    logger.info(f"GPU {device_id}: {free_mem/1e9:.2f}/{total_mem/1e9:.2f} GB free/total")


def clear_gpu_memory():
    """Clear GPU memory pools"""
    if HAS_RMM:
        # RMM handles memory automatically
        pass
    else:
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()
        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
    
    # Force garbage collection
    cp._default_memory_pool.free_all_blocks()
    logger.debug("Cleared GPU memory pools")


def load_preprocessed_data(file_path: Union[str, Path]) -> AnnData:
    """
    Load preprocessed data from h5ad file
    
    Parameters:
    -----------
    file_path : str or Path
        Path to h5ad file
        
    Returns:
    --------
    adata : AnnData
        Loaded data
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"Data file not found: {file_path}")
    
    logger.info(f"Loading data from {file_path}")
    adata = sc.read_h5ad(file_path)
    
    # Verify data is preprocessed
    if 'X_pca' not in adata.obsm:
        logger.warning("No PCA found in data. Data may not be fully preprocessed.")
    
    # Log data info
    logger.info(f"Loaded data: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Check for expected metadata
    if adata.obs.empty:
        logger.warning("No cell metadata found")
    else:
        logger.info(f"Cell metadata columns: {list(adata.obs.columns)}")
    
    return adata


def save_results(results: Dict[str, Any], output_dir: Union[str, Path]):
    """
    Save optimization results
    
    Parameters:
    -----------
    results : dict
        Results dictionary to save
    output_dir : str or Path
        Output directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create timestamped filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save main results as JSON
    results_file = output_dir / f"optimization_results_{timestamp}.json"
    
    # Convert numpy types for JSON serialization
    results_serializable = convert_numpy_types(results)
    
    with open(results_file, 'w') as f:
        json.dump(results_serializable, f, indent=2)
    
    logger.info(f"Saved results to {results_file}")
    
    # Save trials dataframe if present
    if 'all_trials' in results and isinstance(results['all_trials'], dict):
        trials_df = pd.DataFrame(results['all_trials'])
        trials_file = output_dir / f"optimization_trials_{timestamp}.csv"
        trials_df.to_csv(trials_file, index=False)
        logger.info(f"Saved trials to {trials_file}")
    
    # Create summary report
    create_summary_report(results, output_dir / f"optimization_summary_{timestamp}.txt")


def convert_numpy_types(obj):
    """Convert numpy types to Python native types for JSON serialization"""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj


def create_summary_report(results: Dict[str, Any], output_file: Path):
    """Create human-readable summary report"""
    with open(output_file, 'w') as f:
        f.write("LNSC Clustering Optimization Summary\n")
        f.write("=" * 50 + "\n\n")
        
        # Best parameters
        f.write("Best Parameters:\n")
        f.write("-" * 20 + "\n")
        if 'best_params' in results:
            for param, value in results['best_params'].items():
                f.write(f"  {param}: {value}\n")
        f.write("\n")
        
        # Best score
        if 'best_score' in results:
            f.write(f"Best Score: {results['best_score']:.4f}\n\n")
        
        # Study statistics
        if 'study_stats' in results:
            f.write("Optimization Statistics:\n")
            f.write("-" * 20 + "\n")
            for stat, value in results['study_stats'].items():
                f.write(f"  {stat}: {value}\n")
        f.write("\n")
        
        # Validation results if present
        if 'validation' in results:
            f.write("Validation Results:\n")
            f.write("-" * 20 + "\n")
            for sample, scores in results['validation'].items():
                f.write(f"\n  {sample}:\n")
                for metric, value in scores.items():
                    f.write(f"    {metric}: {value:.4f}\n")
    
    logger.info(f"Created summary report: {output_file}")


def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """Load configuration from JSON file"""
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Validate required fields
    required_fields = ['data_path', 'output_dir', 'markers', 'optimization']
    missing_fields = [field for field in required_fields if field not in config]
    
    if missing_fields:
        raise ValueError(f"Missing required configuration fields: {missing_fields}")
    
    return config


def get_gpu_memory_info(device_id: int = 0) -> Dict[str, float]:
    """Get current GPU memory usage information"""
    device = cp.cuda.Device(device_id)
    free_mem, total_mem = device.mem_info
    used_mem = total_mem - free_mem
    
    info = {
        'total_gb': total_mem / 1e9,
        'used_gb': used_mem / 1e9,
        'free_gb': free_mem / 1e9,
        'used_percent': (used_mem / total_mem) * 100
    }
    
    if HAS_RMM:
        # Add RMM pool info if available
        try:
            import rmm
            mr = rmm.mr.get_current_device_resource()
            if hasattr(mr, 'get_info'):
                pool_info = mr.get_info()
                info['pool_size_gb'] = pool_info.get('pool_size', 0) / 1e9
                info['pool_used_gb'] = pool_info.get('used_size', 0) / 1e9
        except Exception:
            pass
    
    return info


def monitor_gpu_memory(func):
    """Decorator to monitor GPU memory usage during function execution"""
    def wrapper(*args, **kwargs):
        # Before execution
        before_info = get_gpu_memory_info()
        logger.debug(f"GPU memory before {func.__name__}: {before_info['used_gb']:.2f} GB used")
        
        try:
            # Execute function
            result = func(*args, **kwargs)
            
            # After execution
            after_info = get_gpu_memory_info()
            memory_delta = after_info['used_gb'] - before_info['used_gb']
            logger.debug(f"GPU memory after {func.__name__}: {after_info['used_gb']:.2f} GB used "
                        f"(delta: {memory_delta:+.2f} GB)")
            
            return result
            
        except cp.cuda.memory.OutOfMemoryError as e:
            logger.error(f"GPU out of memory in {func.__name__}")
            logger.error(f"Memory state: {get_gpu_memory_info()}")
            raise e
    
    return wrapper


def validate_markers(adata: AnnData, markers: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """
    Validate that markers exist in the data
    
    Parameters:
    -----------
    adata : AnnData
        Data to check
    markers : dict
        Dictionary of marker lists
        
    Returns:
    --------
    valid_markers : dict
        Dictionary with only valid markers
    """
    valid_markers = {}
    
    for category, marker_list in markers.items():
        valid_list = []
        
        for marker in marker_list:
            if marker in adata.var_names:
                valid_list.append(marker)
            else:
                # Try case variations
                marker_lower = marker.lower()
                marker_upper = marker.upper()
                marker_title = marker.title()
                
                found = False
                for variant in [marker_lower, marker_upper, marker_title]:
                    if variant in adata.var_names:
                        valid_list.append(variant)
                        logger.warning(f"Found {marker} as {variant}")
                        found = True
                        break
                
                if not found:
                    logger.warning(f"Marker {marker} not found in data")
        
        if valid_list:
            valid_markers[category] = valid_list
        else:
            logger.warning(f"No valid markers found for category {category}")
    
    return valid_markers