# Nextflow Modules

This directory contains modular Nextflow processes and subworkflows for the bioinformatics pipeline.

## Structure

```
modules/
├── local/           # Project-specific modules
│   ├── fastqc/
│   ├── trimming/
│   ├── alignment/
│   └── multiqc/
└── nf-core/         # Modules from nf-core/modules
    └── modules/
```

## Usage

Modules follow the nf-core module standards and can be imported into workflows:

```groovy
include { FASTQC } from './modules/local/fastqc/main'
include { STAR_ALIGN } from './modules/nf-core/modules/star/align/main'
```

## Development

To create a new module:

1. Create a directory with the tool name
2. Add `main.nf` with the process definition
3. Add `meta.yml` with module metadata
4. Add test files in `tests/`

For more information, see the [nf-core modules documentation](https://nf-co.re/docs/contributing/modules).

## Available Modules

- [ ] FASTQC - Quality control of raw reads
- [ ] TRIMMOMATIC - Read trimming and filtering
- [ ] STAR - RNA-seq alignment
- [ ] SALMON - Transcript quantification
- [ ] MULTIQC - Aggregate QC reports
- [ ] SAMTOOLS - SAM/BAM file processing
- [ ] BEDTOOLS - Genomic interval operations

*Note: Modules are placeholders and need to be implemented.*

