#!/usr/bin/env nextflow
/*
========================================================================================
    Bioinformatics Analysis Pipeline
========================================================================================
    Github : https://github.com/bioinfo-project
    
    Main workflow entry point
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow BIOINFO_PIPELINE {
    
    take:
    input_ch
    
    main:
    
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )
    
    //
    // WORKFLOW: Run main analysis
    //
    if (!params.skip_analysis) {
        RUN_ANALYSIS (
            PIPELINE_INITIALISATION.out.samplesheet
        )
    }
    
    emit:
    multiqc_report = RUN_ANALYSIS.out.multiqc_report ?: []
    versions       = RUN_ANALYSIS.out.versions       ?: []
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    main:
    
    //
    // WORKFLOW: Run pipeline
    //
    BIOINFO_PIPELINE (
        Channel.empty()
    )
    
}

/*
========================================================================================
    WORKFLOW: Pipeline initialisation
========================================================================================
*/

workflow PIPELINE_INITIALISATION {
    
    take:
    version
    help
    validate_params
    monochrome_logs
    args
    outdir
    input
    
    main:
    
    // Print version and exit
    if (params.version) {
        log.info "${workflow.manifest.name} v${workflow.manifest.version}"
        System.exit(0)
    }
    
    // Print help and exit
    if (params.help) {
        log.info help_message()
        System.exit(0)
    }
    
    // Validate input parameters
    if (validate_params) {
        validateParameters()
    }
    
    // Check input path parameters
    def checkPathParamList = [
        params.input,
        params.genome,
        params.gtf
    ]
    for (param in checkPathParamList) {
        if (param) {
            file(param, checkIfExists: true)
        }
    }
    
    // Check mandatory parameters
    if (params.input) {
        ch_input = file(params.input)
    } else {
        log.error "Please provide an input samplesheet with the '--input' parameter."
        System.exit(1)
    }
    
    emit:
    samplesheet = ch_input
}

/*
========================================================================================
    WORKFLOW: Main analysis workflow
========================================================================================
*/

workflow RUN_ANALYSIS {
    
    take:
    samplesheet
    
    main:
    
    ch_versions = Channel.empty()
    
    // Parse samplesheet
    ch_samples = samplesheet
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.single_end = row.single_end.toBoolean()
            
            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = meta.single_end ? null : file(row.fastq_2, checkIfExists: true)
            
            return [meta, [fastq_1, fastq_2]].flatten().findAll { it != null }
        }
    
    //
    // MODULE: Quality control
    //
    if (!params.skip_qc) {
        FASTQC (
            ch_samples
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions)
    }
    
    //
    // MODULE: Read trimming
    //
    if (!params.skip_trimming) {
        TRIMMING (
            ch_samples
        )
        ch_trimmed = TRIMMING.out.reads
        ch_versions = ch_versions.mix(TRIMMING.out.versions)
    } else {
        ch_trimmed = ch_samples
    }
    
    //
    // MODULE: Alignment
    //
    if (!params.skip_alignment && params.genome) {
        ALIGNMENT (
            ch_trimmed,
            params.genome
        )
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
    }
    
    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = Channel.empty()
    if (!params.skip_qc) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    }
    
    MULTIQC (
        ch_multiqc_files.collect()
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)
    
    emit:
    multiqc_report = MULTIQC.out.report
    versions       = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def help_message() {
    log.info"""
    Usage:
    
    nextflow run main.nf --input samplesheet.csv --genome genome.fasta [options]
    
    Required arguments:
      --input [file]                  Path to comma-separated file containing information about the samples
      --genome [file]                 Path to FASTA genome file
    
    Optional arguments:
      --outdir [path]                 The output directory where the results will be saved (default: results)
      --gtf [file]                    Path to GTF annotation file
      --skip_qc                       Skip quality control steps
      --skip_trimming                 Skip read trimming
      --skip_alignment                Skip read alignment
      
    Other options:
      --max_memory [memory]           Maximum memory that can be requested for any single job (default: 8.GB)
      --max_cpus [int]                Maximum CPUs that can be requested for any single job (default: 4)
      --max_time [time]               Maximum time that can be requested for any single job (default: 24.h)
      
    """.stripIndent()
}

def validateParameters() {
    // Implementation for parameter validation
    return true
}

/*
========================================================================================
    THE END
========================================================================================
*/

