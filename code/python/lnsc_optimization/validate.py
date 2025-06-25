#!/usr/bin/env python3
"""
Validation script for applying optimized parameters to new samples
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import scanpy as sc

from gpu_clustering import GPUClusterer
from marker_scoring import MarkerScorer
from utils import setup_gpu_memory, load_config, save_results

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_samples(
    params_file: str,
    sample_files: List[str],
    config_file: str,
    output_file: str
):
    """
    Validate optimized parameters on new samples
    
    Parameters:
    -----------
    params_file : str
        JSON file with best parameters
    sample_files : list
        List of h5ad files to validate
    config_file : str
        Configuration file path
    output_file : str
        Output JSON file for results
    """
    # Load configuration
    config = load_config(config_file)
    
    # Setup GPU
    setup_gpu_memory(pool_size_gb=config.get('gpu_memory_gb', 16))
    
    # Load parameters
    with open(params_file, 'r') as f:
        best_params = json.load(f)
    
    logger.info(f"Loaded parameters: {best_params}")
    
    # Initialize modules
    clusterer = GPUClusterer()
    scorer = MarkerScorer(config['markers'])
    
    # Validate each sample
    validation_results = {}
    
    for sample_file in sample_files:
        sample_name = Path(sample_file).stem
        logger.info(f"Validating {sample_name}")
        
        try:
            # Load data
            adata = sc.read_h5ad(sample_file)
            
            # Apply clustering with best parameters
            adata = clusterer.cluster(adata, **best_params)
            
            # Calculate scores
            scores = scorer.calculate_scores(adata)
            
            # Store results
            validation_results[sample_name] = {
                'scores': scores,
                'n_clusters': adata.obs['leiden'].nunique(),
                'n_cells': adata.n_obs,
                'parameters_used': best_params
            }
            
            # Save clustered data
            output_dir = Path(output_file).parent
            output_h5ad = output_dir / f"{sample_name}_validated.h5ad"
            adata.write_h5ad(output_h5ad)
            logger.info(f"Saved clustered data to {output_h5ad}")
            
            # Generate marker report
            report = scorer.generate_marker_report(adata)
            report_file = output_dir / f"{sample_name}_marker_report.csv"
            report.to_csv(report_file, index=False)
            
        except Exception as e:
            logger.error(f"Validation failed for {sample_name}: {str(e)}")
            validation_results[sample_name] = {
                'error': str(e),
                'status': 'failed'
            }
    
    # Save validation results
    with open(output_file, 'w') as f:
        json.dump(validation_results, f, indent=2)
    
    logger.info(f"Validation complete. Results saved to {output_file}")
    
    # Print summary
    print("\nValidation Summary:")
    print("-" * 50)
    for sample, results in validation_results.items():
        print(f"\n{sample}:")
        if 'error' in results:
            print(f"  Status: FAILED - {results['error']}")
        else:
            print(f"  Clusters: {results['n_clusters']}")
            print(f"  Overall Score: {results['scores'].get('overall_score', 0):.3f}")


def main():
    parser = argparse.ArgumentParser(description='Validate optimized parameters on new samples')
    parser.add_argument('--params', required=True, help='JSON file with best parameters')
    parser.add_argument('--samples', nargs='+', required=True, help='Sample h5ad files to validate')
    parser.add_argument('--config', required=True, help='Configuration file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    
    args = parser.parse_args()
    
    validate_samples(
        args.params,
        args.samples,
        args.config,
        args.output
    )


if __name__ == '__main__':
    main()