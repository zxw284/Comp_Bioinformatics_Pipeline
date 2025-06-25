#!/bin/bash
#
# Run LNSC Bayesian Hyperparameter Optimization
# GPU-accelerated clustering optimization for single-cell RNA-seq
#

set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Default values
CONFIG_FILE="${PROJECT_DIR}/config/lnsc_optimization.json"
N_TRIALS=100
N_JOBS=1
PROFILE="gpu"
WORKFLOW="standalone"
DRY_RUN=false
VERBOSE=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check GPU availability
check_gpu() {
    if ! command -v nvidia-smi &> /dev/null; then
        print_error "nvidia-smi not found. Please ensure NVIDIA drivers are installed."
        return 1
    fi
    
    if ! nvidia-smi &> /dev/null; then
        print_error "Failed to communicate with NVIDIA driver."
        return 1
    fi
    
    print_info "GPU check passed"
    nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader
    return 0
}

# Function to check Python dependencies
check_dependencies() {
    print_info "Checking Python dependencies..."
    
    local missing_deps=()
    
    # Check required packages
    for pkg in cupy scanpy optuna pandas numpy; do
        if ! python -c "import $pkg" 2>/dev/null; then
            missing_deps+=("$pkg")
        fi
    done
    
    if [ ${#missing_deps[@]} -gt 0 ]; then
        print_error "Missing Python packages: ${missing_deps[*]}"
        print_info "Please install with: pip install ${missing_deps[*]}"
        return 1
    fi
    
    print_info "All dependencies found"
    return 0
}

# Function to setup environment
setup_environment() {
    print_info "Setting up environment..."
    
    # Export CUDA paths if needed
    if [ -d "/usr/local/cuda" ]; then
        export PATH="/usr/local/cuda/bin:$PATH"
        export LD_LIBRARY_PATH="/usr/local/cuda/lib64:$LD_LIBRARY_PATH"
    fi
    
    # Set Python path
    export PYTHONPATH="${PROJECT_DIR}/code/python:$PYTHONPATH"
    
    # Create output directories
    mkdir -p "${PROJECT_DIR}/results/lnsc_optimization"
    mkdir -p "${PROJECT_DIR}/logs"
}

# Function to run standalone optimization
run_standalone() {
    local data_file="$1"
    local extra_args="${2:-}"
    
    print_info "Running standalone optimization..."
    print_info "Data file: $data_file"
    print_info "Config: $CONFIG_FILE"
    print_info "Trials: $N_TRIALS"
    
    if [ "$DRY_RUN" = true ]; then
        print_warning "DRY RUN - would execute:"
        echo "python ${PROJECT_DIR}/code/python/lnsc_optimization/bayesian_optimizer.py \\"
        echo "    --config $CONFIG_FILE \\"
        echo "    --n-trials $N_TRIALS \\"
        echo "    --n-jobs $N_JOBS \\"
        echo "    $extra_args"
        return 0
    fi
    
    # Update config with data path
    local temp_config=$(mktemp)
    python -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
config['data_path'] = '$data_file'
with open('$temp_config', 'w') as f:
    json.dump(config, f, indent=2)
"
    
    # Run optimization
    python "${PROJECT_DIR}/code/python/lnsc_optimization/bayesian_optimizer.py" \
        --config "$temp_config" \
        --n-trials "$N_TRIALS" \
        --n-jobs "$N_JOBS" \
        $extra_args
    
    # Cleanup
    rm -f "$temp_config"
}

# Function to run via Nextflow
run_nextflow() {
    local data_file="$1"
    
    print_info "Running optimization via Nextflow..."
    
    local nf_args=""
    if [ "$PROFILE" = "gpu" ]; then
        nf_args="-profile docker,gpu"
    else
        nf_args="-profile docker"
    fi
    
    if [ "$DRY_RUN" = true ]; then
        print_warning "DRY RUN - would execute:"
        echo "nextflow run ${PROJECT_DIR}/workflows/nextflow/main.nf \\"
        echo "    -entry LNSC_HYPERPARAMETER_OPTIMIZATION \\"
        echo "    --input $data_file \\"
        echo "    --config $CONFIG_FILE \\"
        echo "    $nf_args"
        return 0
    fi
    
    nextflow run "${PROJECT_DIR}/workflows/nextflow/main.nf" \
        -entry LNSC_HYPERPARAMETER_OPTIMIZATION \
        --input "$data_file" \
        --config "$CONFIG_FILE" \
        $nf_args
}

# Function to show usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] <input_h5ad_file>

Run LNSC Bayesian Hyperparameter Optimization for single-cell clustering

Arguments:
    input_h5ad_file    Path to preprocessed h5ad file from Seurat

Options:
    -c, --config FILE      Configuration file (default: config/lnsc_optimization.json)
    -n, --n-trials N       Number of optimization trials (default: 100)
    -j, --n-jobs N         Number of parallel jobs (default: 1)
    -w, --workflow TYPE    Workflow type: standalone or nextflow (default: standalone)
    -p, --profile PROFILE  Execution profile for Nextflow (default: gpu)
    -v, --validate FILES   Additional h5ad files to validate (space-separated)
    --dry-run             Show what would be executed without running
    --verbose             Enable verbose output
    -h, --help            Show this help message

Examples:
    # Run optimization on in vivo data
    $0 data/processed/lnsc_invivo.h5ad

    # Run with custom config and more trials
    $0 -c my_config.json -n 200 data/processed/lnsc_invivo.h5ad

    # Run and validate on in vitro samples
    $0 -v "data/processed/lnsc_2d.h5ad data/processed/lnsc_3d.h5ad" data/processed/lnsc_invivo.h5ad

    # Run via Nextflow pipeline
    $0 -w nextflow -p docker,gpu data/processed/lnsc_invivo.h5ad

EOF
}

# Parse command line arguments
VALIDATION_FILES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -n|--n-trials)
            N_TRIALS="$2"
            shift 2
            ;;
        -j|--n-jobs)
            N_JOBS="$2"
            shift 2
            ;;
        -w|--workflow)
            WORKFLOW="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -v|--validate)
            VALIDATION_FILES="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            print_error "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            # Positional argument (input file)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

# Check if input file was provided
if [ -z "${INPUT_FILE:-}" ]; then
    print_error "No input file provided"
    usage
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    print_error "Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    print_error "Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Main execution
print_info "LNSC Bayesian Hyperparameter Optimization"
print_info "========================================="

# Check GPU
if ! check_gpu; then
    print_error "GPU check failed. This pipeline requires GPU support."
    exit 1
fi

# Setup environment
setup_environment

# Check dependencies for standalone mode
if [ "$WORKFLOW" = "standalone" ]; then
    if ! check_dependencies; then
        print_error "Dependency check failed"
        exit 1
    fi
fi

# Prepare validation arguments
EXTRA_ARGS=""
if [ -n "$VALIDATION_FILES" ]; then
    EXTRA_ARGS="--validate $VALIDATION_FILES"
fi

# Run optimization
case "$WORKFLOW" in
    standalone)
        run_standalone "$INPUT_FILE" "$EXTRA_ARGS"
        ;;
    nextflow)
        run_nextflow "$INPUT_FILE"
        ;;
    *)
        print_error "Unknown workflow type: $WORKFLOW"
        usage
        exit 1
        ;;
esac

print_info "Optimization complete!"
print_info "Results saved to: ${PROJECT_DIR}/results/lnsc_optimization/"