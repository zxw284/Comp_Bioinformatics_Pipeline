#!/usr/bin/env nextflow
/*
 * LNSC Bayesian Hyperparameter Optimization Module
 * GPU-accelerated clustering optimization
 */

process LNSC_OPTIMIZATION {
    tag "$meta.id"
    label 'process_gpu'
    label 'process_high_memory'
    
    container 'nvcr.io/nvidia/rapids/rapids:22.12-cuda11.5-runtime-ubuntu20.04-py3.9'
    
    input:
    tuple val(meta), path(h5ad)
    path config
    
    output:
    tuple val(meta), path("*_results.json"), emit: results
    tuple val(meta), path("*_clustered.h5ad"), emit: clustered
    path "*.csv", emit: trials
    path "*.txt", emit: report
    path "logs/*.log", emit: logs
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_trials = task.ext.n_trials ?: 100
    """
    # Setup environment
    export CUDA_VISIBLE_DEVICES=0
    mkdir -p logs
    
    # Run optimization
    python ${projectDir}/code/python/lnsc_optimization/bayesian_optimizer.py \\
        --config ${config} \\
        --n-trials ${n_trials} \\
        ${args} \\
        2>&1 | tee logs/${prefix}_optimization.log
    
    # Rename outputs
    mv results/lnsc_optimization/optimization_results_*.json ${prefix}_results.json
    mv results/lnsc_optimization/*_clustered.h5ad ${prefix}_clustered.h5ad
    mv results/lnsc_optimization/optimization_trials_*.csv ${prefix}_trials.csv
    mv results/lnsc_optimization/optimization_summary_*.txt ${prefix}_summary.txt
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results.json
    touch ${prefix}_clustered.h5ad
    touch ${prefix}_trials.csv
    touch ${prefix}_summary.txt
    mkdir -p logs
    touch logs/${prefix}_optimization.log
    """
}

process LNSC_VALIDATE {
    tag "$meta.id"
    label 'process_gpu'
    label 'process_medium'
    
    container 'nvcr.io/nvidia/rapids/rapids:22.12-cuda11.5-runtime-ubuntu20.04-py3.9'
    
    input:
    tuple val(meta), path(results_json)
    path validation_samples
    path config
    
    output:
    tuple val(meta), path("*_validation.json"), emit: validation
    path "*.h5ad", emit: validated_samples
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract best parameters from results
    python -c "
import json
with open('${results_json}', 'r') as f:
    results = json.load(f)
    best_params = results['best_params']
    with open('best_params.json', 'w') as out:
        json.dump(best_params, out)
    "
    
    # Validate on each sample
    python ${projectDir}/code/python/lnsc_optimization/validate.py \\
        --params best_params.json \\
        --samples ${validation_samples} \\
        --config ${config} \\
        --output ${prefix}_validation.json
    """
}

workflow LNSC_HYPERPARAMETER_OPTIMIZATION {
    take:
    ch_input        // channel: [ val(meta), path(h5ad) ]
    ch_config       // channel: path(config)
    ch_validation   // channel: path(validation_samples) (optional)
    
    main:
    ch_versions = Channel.empty()
    
    // Run optimization
    LNSC_OPTIMIZATION(
        ch_input,
        ch_config
    )
    
    // Validate if samples provided
    if (ch_validation) {
        LNSC_VALIDATE(
            LNSC_OPTIMIZATION.out.results,
            ch_validation,
            ch_config
        )
        ch_validated = LNSC_VALIDATE.out.validated_samples
    } else {
        ch_validated = Channel.empty()
    }
    
    emit:
    results     = LNSC_OPTIMIZATION.out.results
    clustered   = LNSC_OPTIMIZATION.out.clustered
    trials      = LNSC_OPTIMIZATION.out.trials
    report      = LNSC_OPTIMIZATION.out.report
    validated   = ch_validated
    versions    = ch_versions
}