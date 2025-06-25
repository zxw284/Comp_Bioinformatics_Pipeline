#!/usr/bin/env python3
"""
Marker-based scoring functions for biological validation of clusters
"""

import logging
from typing import Dict, List, Tuple, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

logger = logging.getLogger(__name__)


class MarkerScorer:
    """Score clustering results based on biological markers"""
    
    def __init__(self, marker_config: Dict):
        """
        Initialize with marker configuration
        
        Parameters:
        -----------
        marker_config : dict
            Configuration containing tier1_markers, gradient_markers, etc.
        """
        self.tier1_markers = marker_config.get('tier1_markers', {})
        self.gradient_markers = marker_config.get('gradient_markers', {})
        self.thresholds = marker_config.get('thresholds', {
            'min_fold_change': 2.0,
            'max_pvalue': 0.05,
            'min_pct_expressing': 0.25
        })
        
    def calculate_scores(self, adata: AnnData) -> Dict[str, float]:
        """
        Calculate all scoring metrics
        
        Parameters:
        -----------
        adata : AnnData
            Clustered data with 'leiden' in obs
            
        Returns:
        --------
        scores : dict
            Dictionary of score names to values
        """
        if 'leiden' not in adata.obs:
            raise ValueError("No clustering found. Need 'leiden' in adata.obs")
        
        scores = {}
        
        # Tier 1 marker specificity
        tier1_score = self.score_tier1_markers(adata)
        scores['tier1_specificity'] = tier1_score
        
        # Gradient marker separation
        gradient_score = self.score_gradient_markers(adata)
        scores['gradient_separation'] = gradient_score
        
        # Marker coverage (what fraction of expected markers are well-separated)
        coverage_score = self.score_marker_coverage(adata)
        scores['marker_coverage'] = coverage_score
        
        # Rare population detection
        rare_score = self.score_rare_populations(adata)
        scores['rare_population_score'] = rare_score
        
        # Cluster stability metrics
        stability_score = self.score_cluster_stability(adata)
        scores['stability'] = stability_score
        
        # Overall biological coherence
        scores['overall_score'] = self._calculate_overall_score(scores)
        
        logger.info(f"Calculated scores: {scores}")
        
        return scores
    
    def score_tier1_markers(self, adata: AnnData) -> float:
        """
        Score based on tier 1 marker specificity
        
        High score if each cell type's markers are highly specific to one cluster
        """
        specificities = []
        
        for cell_type, markers in self.tier1_markers.items():
            # Get valid markers (present in data)
            valid_markers = [m for m in markers if m in adata.var_names]
            
            if not valid_markers:
                logger.warning(f"No valid markers found for {cell_type}")
                continue
            
            # Calculate mean expression per cluster
            cluster_expr = self._get_cluster_expression(adata, valid_markers)
            
            # Calculate specificity score
            specificity = self._calculate_marker_specificity(cluster_expr)
            specificities.append(specificity)
            
            logger.debug(f"{cell_type} specificity: {specificity:.3f}")
        
        # Return mean specificity across all cell types
        return np.mean(specificities) if specificities else 0.0
    
    def score_gradient_markers(self, adata: AnnData) -> float:
        """
        Score based on gradient marker separation
        
        High score if gradient markers show clear high/low populations
        """
        gradient_scores = []
        
        for marker_group, marker_info in self.gradient_markers.items():
            if marker_group == 'TRC_subtypes':
                # Special handling for CCL19 high/low
                score = self._score_ccl19_gradient(adata, marker_info)
                if score is not None:
                    gradient_scores.append(score)
            else:
                # Handle other gradient markers
                markers = marker_info if isinstance(marker_info, list) else [marker_info]
                for marker in markers:
                    if marker in adata.var_names:
                        score = self._score_single_gradient(adata, marker)
                        if score is not None:
                            gradient_scores.append(score)
        
        return np.mean(gradient_scores) if gradient_scores else 0.0
    
    def score_marker_coverage(self, adata: AnnData) -> float:
        """
        Score based on how many expected markers are well-separated
        """
        total_markers = 0
        well_separated = 0
        
        # Check all tier 1 markers
        for cell_type, markers in self.tier1_markers.items():
            for marker in markers:
                if marker in adata.var_names:
                    total_markers += 1
                    if self._is_marker_well_separated(adata, marker):
                        well_separated += 1
        
        # Check gradient markers
        for marker_group, marker_info in self.gradient_markers.items():
            if isinstance(marker_info, dict):
                # Complex gradient marker
                for marker in marker_info.keys():
                    if marker in adata.var_names:
                        total_markers += 1
                        if self._is_marker_well_separated(adata, marker):
                            well_separated += 1
            else:
                # Simple list of markers
                markers = marker_info if isinstance(marker_info, list) else [marker_info]
                for marker in markers:
                    if marker in adata.var_names:
                        total_markers += 1
                        if self._is_marker_well_separated(adata, marker):
                            well_separated += 1
        
        return well_separated / total_markers if total_markers > 0 else 0.0
    
    def score_rare_populations(self, adata: AnnData) -> float:
        """
        Score based on detection of rare populations (e.g., Nr4a1+ cells)
        """
        rare_markers = ['Nr4a1', 'Cxcl9', 'Cxcl10']
        scores = []
        
        for marker in rare_markers:
            if marker not in adata.var_names:
                continue
            
            # Check if any cluster is enriched for this marker
            cluster_expr = self._get_cluster_expression(adata, [marker])
            
            # Look for clusters with high expression
            max_expr = cluster_expr[marker].max()
            mean_expr = cluster_expr[marker].mean()
            
            if max_expr > 0 and mean_expr > 0:
                # Score based on fold enrichment of highest cluster
                enrichment = max_expr / mean_expr
                scores.append(min(enrichment / 3.0, 1.0))  # Normalize to 0-1
        
        return np.mean(scores) if scores else 0.0
    
    def score_cluster_stability(self, adata: AnnData) -> float:
        """
        Score cluster stability using bootstrapping
        """
        # For efficiency, we'll use a simple metric based on cluster sizes
        # and compactness in PCA space
        
        cluster_sizes = adata.obs['leiden'].value_counts()
        n_clusters = len(cluster_sizes)
        
        # Penalize too many or too few clusters
        if n_clusters < 5:
            size_penalty = 0.5
        elif n_clusters > 20:
            size_penalty = 0.5
        else:
            size_penalty = 1.0
        
        # Check cluster size distribution (prefer balanced clusters)
        size_cv = cluster_sizes.std() / cluster_sizes.mean()
        size_score = 1.0 / (1.0 + size_cv)
        
        # Check cluster compactness in PCA space
        compactness_scores = []
        X_pca = adata.obsm['X_pca']
        
        for cluster in adata.obs['leiden'].unique():
            mask = adata.obs['leiden'] == cluster
            cluster_pca = X_pca[mask]
            
            if len(cluster_pca) > 10:
                # Calculate within-cluster variance
                cluster_var = np.var(cluster_pca, axis=0).sum()
                compactness_scores.append(1.0 / (1.0 + cluster_var))
        
        compactness = np.mean(compactness_scores) if compactness_scores else 0.5
        
        return size_penalty * size_score * compactness
    
    def _get_cluster_expression(self, adata: AnnData, markers: List[str]) -> pd.DataFrame:
        """Get mean expression of markers per cluster"""
        expr_data = {}
        
        for marker in markers:
            if marker not in adata.var_names:
                continue
            
            marker_idx = adata.var_names.get_loc(marker)
            cluster_means = []
            
            for cluster in sorted(adata.obs['leiden'].unique()):
                mask = adata.obs['leiden'] == cluster
                if hasattr(adata.X, 'toarray'):
                    cluster_expr = adata.X[mask, marker_idx].toarray().flatten()
                else:
                    cluster_expr = adata.X[mask, marker_idx].flatten()
                cluster_means.append(np.mean(cluster_expr))
            
            expr_data[marker] = cluster_means
        
        return pd.DataFrame(expr_data, index=sorted(adata.obs['leiden'].unique()))
    
    def _calculate_marker_specificity(self, cluster_expr: pd.DataFrame) -> float:
        """Calculate how specific markers are to clusters"""
        if cluster_expr.empty:
            return 0.0
        
        specificities = []
        
        for marker in cluster_expr.columns:
            expr_values = cluster_expr[marker].values
            
            if expr_values.max() == 0:
                continue
            
            # Calculate specificity as ratio of max to second max
            sorted_expr = np.sort(expr_values)[::-1]
            if len(sorted_expr) > 1 and sorted_expr[1] > 0:
                specificity = sorted_expr[0] / sorted_expr[1]
            else:
                specificity = sorted_expr[0] / (np.mean(expr_values) + 1e-6)
            
            specificities.append(min(specificity / self.thresholds['min_fold_change'], 1.0))
        
        return np.mean(specificities) if specificities else 0.0
    
    def _score_ccl19_gradient(self, adata: AnnData, marker_info: Dict) -> Optional[float]:
        """Score CCL19 high/low gradient separation"""
        if 'Ccl19' not in adata.var_names:
            return None
        
        ccl19_info = marker_info.get('Ccl19', {})
        high_threshold = ccl19_info.get('high', 2.5)
        low_range = ccl19_info.get('low', (1.0, 2.5))
        
        # Get Ccl19 expression per cluster
        cluster_expr = self._get_cluster_expression(adata, ['Ccl19'])
        ccl19_values = cluster_expr['Ccl19'].values
        
        # Check if we have clear high and low populations
        high_clusters = np.sum(ccl19_values > high_threshold)
        low_clusters = np.sum((ccl19_values >= low_range[0]) & (ccl19_values <= low_range[1]))
        
        if high_clusters > 0 and low_clusters > 0:
            # Good separation
            separation = (ccl19_values.max() - ccl19_values.min()) / (ccl19_values.mean() + 1e-6)
            return min(separation / 2.0, 1.0)
        else:
            return 0.3  # Partial credit if some separation exists
    
    def _score_single_gradient(self, adata: AnnData, marker: str) -> Optional[float]:
        """Score a single gradient marker"""
        cluster_expr = self._get_cluster_expression(adata, [marker])
        if marker not in cluster_expr.columns:
            return None
        
        values = cluster_expr[marker].values
        
        # Check for gradient (variance in expression)
        if values.std() > 0:
            cv = values.std() / (values.mean() + 1e-6)
            return min(cv, 1.0)
        else:
            return 0.0
    
    def _is_marker_well_separated(self, adata: AnnData, marker: str) -> bool:
        """Check if a marker shows good separation between clusters"""
        cluster_expr = self._get_cluster_expression(adata, [marker])
        if marker not in cluster_expr.columns:
            return False
        
        values = cluster_expr[marker].values
        
        # Check if max expression is significantly higher than mean
        if values.max() > 0:
            fold_change = values.max() / (values.mean() + 1e-6)
            return fold_change >= self.thresholds['min_fold_change']
        
        return False
    
    def _calculate_overall_score(self, scores: Dict[str, float]) -> float:
        """Calculate weighted overall score"""
        # Define weights for different components
        weights = {
            'tier1_specificity': 0.35,
            'gradient_separation': 0.25,
            'marker_coverage': 0.20,
            'rare_population_score': 0.10,
            'stability': 0.10
        }
        
        overall = 0.0
        total_weight = 0.0
        
        for key, weight in weights.items():
            if key in scores and scores[key] is not None:
                overall += weight * scores[key]
                total_weight += weight
        
        return overall / total_weight if total_weight > 0 else 0.0
    
    def generate_marker_report(self, adata: AnnData) -> pd.DataFrame:
        """Generate detailed marker expression report"""
        report_data = []
        
        # Analyze each marker
        all_markers = []
        for markers in self.tier1_markers.values():
            all_markers.extend(markers)
        
        for marker in set(all_markers):
            if marker not in adata.var_names:
                continue
            
            # Get expression stats
            marker_idx = adata.var_names.get_loc(marker)
            
            for cluster in sorted(adata.obs['leiden'].unique()):
                mask = adata.obs['leiden'] == cluster
                if hasattr(adata.X, 'toarray'):
                    cluster_expr = adata.X[mask, marker_idx].toarray().flatten()
                else:
                    cluster_expr = adata.X[mask, marker_idx].flatten()
                
                report_data.append({
                    'marker': marker,
                    'cluster': cluster,
                    'mean_expression': np.mean(cluster_expr),
                    'pct_expressing': np.mean(cluster_expr > 0) * 100,
                    'n_cells': mask.sum()
                })
        
        return pd.DataFrame(report_data)