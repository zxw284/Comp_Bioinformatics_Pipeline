#!/bin/bash
# Bioinformatics Analysis Runner Script
# This script runs bioinformatics analysis workflows

set -euo pipefail

# Script metadata
SCRIPT_NAME="run_analysis.sh"
SCRIPT_VERSION="1.0.0"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Logging functions
log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }

# Default values
WORKFLOW="nextflow"
PROFILE="docker"
INPUT="data/samples.tsv"
OUTDIR="results"
DRY_RUN=false
VERBOSE=false
FORCE=false

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Run bioinformatics analysis workflows

OPTIONS:
    -h, --help              Show this help message
    -w, --workflow TYPE     Workflow to run (nextflow, snakemake) [default: nextflow]
    -p, --profile PROFILE   Execution profile (docker, singularity, conda, hpc) [default: docker]
    -i, --input FILE        Input samplesheet [default: data/samples.tsv]
    -o, --outdir DIR        Output directory [default: results]
    -n, --dry-run          Perform dry run without execution
    -v, --verbose          Enable verbose output
    -f, --force            Force overwrite existing results
    --cores N              Number of CPU cores to use [default: 4]
    --memory SIZE          Maximum memory to use [default: 8GB]
    --config FILE          Custom configuration file

EXAMPLES:
    $0                      # Run default Nextflow pipeline
    $0 -w snakemake         # Run Snakemake workflow
    $0 -p conda             # Use conda environments
    $0 --dry-run            # Show what would be executed

EOF
}

# Parse command line arguments
CORES=4
MEMORY="8GB"
CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -w|--workflow)
            WORKFLOW="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --memory)
            MEMORY="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Set project root
PROJ_ROOT=${PROJ_ROOT:-$(pwd)}
cd "$PROJ_ROOT"

# Validate inputs
validate_inputs() {
    log_info "Validating inputs..."
    
    # Check input file
    if [[ ! -f "$INPUT" ]]; then
        log_error "Input file not found: $INPUT"
        exit 1
    fi
    
    # Check workflow type
    if [[ ! "$WORKFLOW" =~ ^(nextflow|snakemake)$ ]]; then
        log_error "Invalid workflow type: $WORKFLOW"
        log_error "Supported workflows: nextflow, snakemake"
        exit 1
    fi
    
    # Check if workflow files exist
    case $WORKFLOW in
        nextflow)
            if [[ ! -f "workflows/nextflow/main.nf" ]]; then
                log_error "Nextflow main.nf not found"
                exit 1
            fi
            ;;
        snakemake)
            if [[ ! -f "workflows/snakemake/Snakefile" ]]; then
                log_error "Snakemake Snakefile not found"
                exit 1
            fi
            ;;
    esac
    
    log_success "Input validation passed"
}

# Check dependencies
check_dependencies() {
    log_info "Checking workflow dependencies..."
    
    case $WORKFLOW in
        nextflow)
            if ! command -v nextflow >/dev/null 2>&1; then
                log_error "Nextflow not found. Please install Nextflow first."
                exit 1
            fi
            ;;
        snakemake)
            if ! command -v snakemake >/dev/null 2>&1; then
                log_error "Snakemake not found. Please install Snakemake first."
                exit 1
            fi
            ;;
    esac
    
    # Check profile dependencies
    case $PROFILE in
        docker)
            if ! command -v docker >/dev/null 2>&1; then
                log_warning "Docker not found. Consider using conda profile."
            fi
            ;;
        singularity)
            if ! command -v singularity >/dev/null 2>&1; then
                log_warning "Singularity not found. Consider using conda profile."
            fi
            ;;
        conda)
            if ! command -v conda >/dev/null 2>&1; then
                log_error "Conda not found. Please install conda first."
                exit 1
            fi
            ;;
    esac
    
    log_success "Dependencies check passed"
}

# Prepare output directory
prepare_output() {
    log_info "Preparing output directory: $OUTDIR"
    
    if [[ -d "$OUTDIR" ]] && [[ "$FORCE" = false ]]; then
        log_warning "Output directory already exists: $OUTDIR"
        read -p "Continue and potentially overwrite results? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "Analysis cancelled"
            exit 0
        fi
    fi
    
    mkdir -p "$OUTDIR"
    log_success "Output directory prepared"
}

# Run Nextflow workflow
run_nextflow() {
    log_info "Running Nextflow workflow..."
    
    local cmd=("nextflow" "run" "workflows/nextflow/main.nf")
    
    # Add parameters
    cmd+=("--input" "$INPUT")
    cmd+=("--outdir" "$OUTDIR")
    
    # Add profile
    cmd+=("-profile" "$PROFILE")
    
    # Add resources
    cmd+=("--max_cpus" "$CORES")
    cmd+=("--max_memory" "$MEMORY")
    
    # Add configuration
    if [[ -n "$CONFIG_FILE" ]]; then
        cmd+=("-c" "$CONFIG_FILE")
    fi
    
    # Add optional flags
    if [[ "$VERBOSE" = true ]]; then
        cmd+=("-ansi-log" "false")
    fi
    
    # Dry run
    if [[ "$DRY_RUN" = true ]]; then
        log_info "Dry run mode - showing command that would be executed:"
        echo "${cmd[*]}"
        return 0
    fi
    
    # Execute command
    log_info "Executing: ${cmd[*]}"
    "${cmd[@]}"
}

# Run Snakemake workflow
run_snakemake() {
    log_info "Running Snakemake workflow..."
    
    local cmd=("snakemake")
    
    # Add Snakefile
    cmd+=("--snakefile" "workflows/snakemake/Snakefile")
    
    # Add cores
    cmd+=("--cores" "$CORES")
    
    # Add configuration
    cmd+=("--configfile" "config/project.yml")
    
    # Add profile-specific options
    case $PROFILE in
        conda)
            cmd+=("--use-conda")
            ;;
        singularity)
            cmd+=("--use-singularity")
            ;;
    esac
    
    # Add other options
    cmd+=("--printshellcmds")
    cmd+=("--reason")
    
    if [[ "$VERBOSE" = true ]]; then
        cmd+=("--verbose")
    fi
    
    if [[ "$FORCE" = true ]]; then
        cmd+=("--forceall")
    fi
    
    # Dry run
    if [[ "$DRY_RUN" = true ]]; then
        cmd+=("--dry-run")
        log_info "Dry run mode enabled"
    fi
    
    # Execute command
    log_info "Executing: ${cmd[*]}"
    cd workflows/snakemake
    "${cmd[@]}"
    cd "$PROJ_ROOT"
}

# Generate report
generate_report() {
    log_info "Generating analysis report..."
    
    local report_dir="$OUTDIR/reports"
    mkdir -p "$report_dir"
    
    # Create summary report
    cat > "$report_dir/analysis_summary.md" << EOF
# Analysis Summary

**Date:** $(date)
**Workflow:** $WORKFLOW
**Profile:** $PROFILE
**Input:** $INPUT
**Output:** $OUTDIR

## Configuration
- Cores: $CORES
- Memory: $MEMORY
- Force: $FORCE
- Dry Run: $DRY_RUN

## Files Generated

EOF
    
    # List output files
    if [[ -d "$OUTDIR" ]]; then
        echo "\`\`\`" >> "$report_dir/analysis_summary.md"
        find "$OUTDIR" -type f -name "*.html" -o -name "*.pdf" -o -name "*.png" | head -20 >> "$report_dir/analysis_summary.md"
        echo "\`\`\`" >> "$report_dir/analysis_summary.md"
    fi
    
    log_success "Report generated: $report_dir/analysis_summary.md"
}

# Main execution
main() {
    log_info "Starting $WORKFLOW analysis..."
    log_info "Profile: $PROFILE"
    log_info "Input: $INPUT"
    log_info "Output: $OUTDIR"
    
    # Validation
    validate_inputs
    check_dependencies
    
    # Preparation
    if [[ "$DRY_RUN" = false ]]; then
        prepare_output
    fi
    
    # Run workflow
    case $WORKFLOW in
        nextflow)
            run_nextflow
            ;;
        snakemake)
            run_snakemake
            ;;
    esac
    
    # Post-processing
    if [[ "$DRY_RUN" = false ]]; then
        generate_report
        log_success "Analysis completed successfully!"
        log_info "Results available in: $OUTDIR"
    else
        log_success "Dry run completed"
    fi
}

# Error handling
trap 'log_error "Analysis failed at line $LINENO"' ERR

# Execute main function
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi

