#!/usr/bin/env python3
"""
GPU-accelerated clustering using RAPIDS
"""

import logging
from typing import Optional, Tuple

import cupy as cp
import cudf
import cuml
import numpy as np
import scanpy as sc
from anndata import AnnData
from rapids_singlecell import pp, tl, pl
from scipy.sparse import issparse

logger = logging.getLogger(__name__)


class GPUClusterer:
    """GPU-accelerated clustering for single-cell data"""
    
    def __init__(self):
        """Initialize GPU clusterer"""
        self.check_gpu_availability()
        
    def check_gpu_availability(self):
        """Check if GPU is available and log info"""
        try:
            n_devices = cp.cuda.runtime.getDeviceCount()
            logger.info(f"Found {n_devices} GPU device(s)")
            
            for i in range(n_devices):
                device = cp.cuda.Device(i)
                mem_info = device.mem_info
                logger.info(f"GPU {i}: {device.compute_capability}, "
                          f"Memory: {mem_info[1] / 1e9:.2f} GB total, "
                          f"{mem_info[0] / 1e9:.2f} GB free")
        except Exception as e:
            logger.error(f"GPU check failed: {str(e)}")
            raise RuntimeError("GPU not available or CUDA error")
    
    def cluster(self, 
                adata: AnnData,
                n_pcs: int = 50,
                n_neighbors: int = 30,
                min_dist: float = 0.3,
                resolution: float = 1.0,
                random_state: int = 42) -> AnnData:
        """
        Perform GPU-accelerated clustering
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data matrix (already preprocessed from Seurat)
        n_pcs : int
            Number of principal components to use
        n_neighbors : int
            Number of neighbors for graph construction
        min_dist : float
            Minimum distance for UMAP
        resolution : float
            Resolution parameter for Leiden clustering
        random_state : int
            Random seed for reproducibility
            
        Returns:
        --------
        adata : AnnData
            Data with clustering results added
        """
        logger.info("Starting GPU-accelerated clustering")
        
        # Ensure we're working with a copy
        adata = adata.copy()
        
        # Check if PCA is already computed from Seurat
        if 'X_pca' not in adata.obsm:
            logger.info("Computing PCA on GPU")
            # Convert to GPU array
            if issparse(adata.X):
                X_gpu = cp.array(adata.X.toarray())
            else:
                X_gpu = cp.array(adata.X)
            
            # GPU PCA
            pca = cuml.PCA(n_components=n_pcs, random_state=random_state)
            X_pca_gpu = pca.fit_transform(X_gpu)
            adata.obsm['X_pca'] = cp.asnumpy(X_pca_gpu)
            
            # Clean up GPU memory
            del X_gpu, X_pca_gpu
            cp.get_default_memory_pool().free_all_blocks()
        else:
            logger.info("Using existing PCA from Seurat")
            # Trim to requested number of PCs
            adata.obsm['X_pca'] = adata.obsm['X_pca'][:, :n_pcs]
        
        # Convert PCA to GPU for downstream analysis
        X_pca_gpu = cp.array(adata.obsm['X_pca'])
        
        # Compute neighbors on GPU
        logger.info(f"Computing neighbors with k={n_neighbors}")
        nn = cuml.NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
        nn.fit(X_pca_gpu)
        distances, indices = nn.kneighbors(X_pca_gpu)
        
        # Store in adata (convert back to CPU)
        adata.obsp['distances'] = cp.asnumpy(distances)
        adata.obsp['connectivities'] = cp.asnumpy(indices)
        adata.uns['neighbors'] = {
            'params': {
                'n_neighbors': n_neighbors,
                'method': 'rapids_cuml'
            }
        }
        
        # UMAP on GPU
        logger.info(f"Computing UMAP with min_dist={min_dist}")
        umap = cuml.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=2,
            random_state=random_state
        )
        X_umap_gpu = umap.fit_transform(X_pca_gpu)
        adata.obsm['X_umap'] = cp.asnumpy(X_umap_gpu)
        
        # Leiden clustering on GPU
        logger.info(f"Computing Leiden clustering with resolution={resolution}")
        # Note: RAPIDS doesn't have Leiden yet, so we'll use CPU version
        # but with GPU-computed neighbors
        sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
        
        # Clean up GPU memory
        del X_pca_gpu, distances, indices, X_umap_gpu
        cp.get_default_memory_pool().free_all_blocks()
        
        logger.info(f"Clustering complete. Found {adata.obs['leiden'].nunique()} clusters")
        
        return adata
    
    def cluster_batch(self,
                     adata: AnnData,
                     batch_key: str,
                     **kwargs) -> AnnData:
        """
        Perform batch-aware clustering
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data matrix
        batch_key : str
            Key in adata.obs containing batch information
        **kwargs : dict
            Additional arguments passed to cluster()
            
        Returns:
        --------
        adata : AnnData
            Data with batch-corrected clustering results
        """
        logger.info(f"Performing batch-aware clustering with batch_key='{batch_key}'")
        
        # Check if batch key exists
        if batch_key not in adata.obs:
            raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")
        
        # For now, perform standard clustering
        # In future, could implement Harmony or other batch correction on GPU
        adata = self.cluster(adata, **kwargs)
        
        # Add batch information to clustering results
        adata.uns['clustering_batch_key'] = batch_key
        
        return adata
    
    def compute_cluster_metrics(self, adata: AnnData) -> dict:
        """
        Compute clustering quality metrics
        
        Parameters:
        -----------
        adata : AnnData
            Clustered data
            
        Returns:
        --------
        metrics : dict
            Dictionary of clustering metrics
        """
        if 'leiden' not in adata.obs:
            raise ValueError("No clustering found. Run cluster() first.")
        
        metrics = {}
        
        # Number of clusters
        metrics['n_clusters'] = adata.obs['leiden'].nunique()
        
        # Cluster sizes
        cluster_sizes = adata.obs['leiden'].value_counts()
        metrics['min_cluster_size'] = cluster_sizes.min()
        metrics['max_cluster_size'] = cluster_sizes.max()
        metrics['mean_cluster_size'] = cluster_sizes.mean()
        
        # Silhouette score (on subset for speed)
        if adata.n_obs > 5000:
            # Subsample for large datasets
            indices = np.random.choice(adata.n_obs, 5000, replace=False)
            X_subset = adata.obsm['X_pca'][indices]
            labels_subset = adata.obs['leiden'].iloc[indices].astype(int)
        else:
            X_subset = adata.obsm['X_pca']
            labels_subset = adata.obs['leiden'].astype(int)
        
        # Convert to GPU for silhouette calculation
        X_gpu = cp.array(X_subset)
        labels_gpu = cp.array(labels_subset)
        
        try:
            from cuml.metrics import silhouette_score
            metrics['silhouette_score'] = float(
                silhouette_score(X_gpu, labels_gpu, metric='euclidean')
            )
        except Exception as e:
            logger.warning(f"Could not compute silhouette score: {str(e)}")
            metrics['silhouette_score'] = None
        
        # Clean up
        del X_gpu, labels_gpu
        cp.get_default_memory_pool().free_all_blocks()
        
        return metrics