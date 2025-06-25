#!/usr/bin/env python3
"""
Bayesian Hyperparameter Optimization for LNSC Clustering
GPU-accelerated optimization using RAPIDS and Optuna
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import cupy as cp
import numpy as np
import optuna
import pandas as pd
import scanpy as sc
from optuna.trial import Trial
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner

from gpu_clustering import GPUClusterer
from marker_scoring import MarkerScorer
from utils import setup_gpu_memory, load_preprocessed_data, save_results


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class LNSCOptimizer:
    """Bayesian optimizer for LNSC clustering parameters"""
    
    def __init__(self, config_path: str):
        """Initialize optimizer with configuration"""
        with open(config_path, 'r') as f:
            self.config = json.load(f)
        
        self.clusterer = GPUClusterer()
        self.scorer = MarkerScorer(self.config['markers'])
        self.best_params = None
        self.best_score = -np.inf
        
        # Setup GPU memory
        setup_gpu_memory(pool_size_gb=self.config.get('gpu_memory_gb', 16))
        
    def objective(self, trial: Trial) -> float:
        """Objective function for Optuna optimization"""
        try:
            # Sample hyperparameters
            params = self._sample_parameters(trial)
            
            # Load data (already preprocessed from Seurat)
            adata = load_preprocessed_data(self.config['data_path'])
            
            # Run clustering with sampled parameters
            adata = self.clusterer.cluster(
                adata,
                n_pcs=params['n_pcs'],
                n_neighbors=params['n_neighbors'],
                min_dist=params['min_dist'],
                resolution=params['resolution'],
                random_state=params.get('random_state', 42)
            )
            
            # Calculate marker-based scores
            scores = self.scorer.calculate_scores(adata)
            
            # Composite score (weighted average)
            composite_score = self._calculate_composite_score(scores)
            
            # Report intermediate values for pruning
            trial.report(composite_score, step=0)
            
            # Clear GPU memory
            cp.get_default_memory_pool().free_all_blocks()
            
            return composite_score
            
        except Exception as e:
            logger.error(f"Trial failed: {str(e)}")
            return -np.inf
    
    def _sample_parameters(self, trial: Trial) -> Dict[str, Any]:
        """Sample hyperparameters for trial"""
        param_ranges = self.config['optimization']['param_ranges']
        
        params = {
            'n_pcs': trial.suggest_int('n_pcs', 
                                      param_ranges['n_pcs']['min'], 
                                      param_ranges['n_pcs']['max']),
            'n_neighbors': trial.suggest_int('n_neighbors', 
                                           param_ranges['n_neighbors']['min'], 
                                           param_ranges['n_neighbors']['max']),
            'min_dist': trial.suggest_float('min_dist', 
                                          param_ranges['min_dist']['min'], 
                                          param_ranges['min_dist']['max'], 
                                          log=True),
            'resolution': trial.suggest_float('resolution', 
                                            param_ranges['resolution']['min'], 
                                            param_ranges['resolution']['max'], 
                                            log=True),
            'random_state': self.config.get('random_state', 42)
        }
        
        return params
    
    def _calculate_composite_score(self, scores: Dict[str, float]) -> float:
        """Calculate weighted composite score from individual metrics"""
        weights = self.config['optimization']['score_weights']
        
        composite = 0.0
        for metric, weight in weights.items():
            if metric in scores:
                composite += weight * scores[metric]
        
        return composite
    
    def optimize(self, n_trials: int = 100, n_jobs: int = 1) -> Dict[str, Any]:
        """Run Bayesian optimization"""
        logger.info(f"Starting optimization with {n_trials} trials")
        
        # Create study
        study = optuna.create_study(
            direction='maximize',
            sampler=TPESampler(seed=self.config.get('random_state', 42)),
            pruner=MedianPruner(n_startup_trials=10, n_warmup_steps=5),
            study_name='lnsc_clustering_optimization'
        )
        
        # Optimize
        study.optimize(
            self.objective,
            n_trials=n_trials,
            n_jobs=n_jobs,
            show_progress_bar=True
        )
        
        # Get best parameters
        self.best_params = study.best_params
        self.best_score = study.best_value
        
        logger.info(f"Best score: {self.best_score:.4f}")
        logger.info(f"Best parameters: {self.best_params}")
        
        # Save results
        results = {
            'best_params': self.best_params,
            'best_score': self.best_score,
            'study_stats': {
                'n_trials': len(study.trials),
                'n_pruned': len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED]),
                'n_complete': len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
            },
            'all_trials': study.trials_dataframe().to_dict()
        }
        
        save_results(results, self.config['output_dir'])
        
        return results
    
    def validate_on_samples(self, sample_paths: List[str]) -> Dict[str, Dict[str, float]]:
        """Validate optimized parameters on additional samples"""
        logger.info("Validating on additional samples")
        
        if self.best_params is None:
            raise ValueError("No optimized parameters found. Run optimize() first.")
        
        validation_results = {}
        
        for sample_path in sample_paths:
            sample_name = Path(sample_path).stem
            logger.info(f"Validating on {sample_name}")
            
            # Load and cluster
            adata = load_preprocessed_data(sample_path)
            adata = self.clusterer.cluster(adata, **self.best_params)
            
            # Score
            scores = self.scorer.calculate_scores(adata)
            validation_results[sample_name] = scores
            
            # Save clustered data
            output_path = Path(self.config['output_dir']) / f"{sample_name}_clustered.h5ad"
            adata.write_h5ad(output_path)
            
        return validation_results


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='LNSC Bayesian Hyperparameter Optimization')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--n-trials', type=int, default=100, help='Number of optimization trials')
    parser.add_argument('--n-jobs', type=int, default=1, help='Number of parallel jobs')
    parser.add_argument('--validate', nargs='+', help='Additional samples for validation')
    
    args = parser.parse_args()
    
    # Initialize optimizer
    optimizer = LNSCOptimizer(args.config)
    
    # Run optimization
    results = optimizer.optimize(n_trials=args.n_trials, n_jobs=args.n_jobs)
    
    # Validate on additional samples if provided
    if args.validate:
        validation_results = optimizer.validate_on_samples(args.validate)
        results['validation'] = validation_results
        
        # Save updated results
        save_results(results, optimizer.config['output_dir'])
    
    logger.info("Optimization complete!")


if __name__ == '__main__':
    main()