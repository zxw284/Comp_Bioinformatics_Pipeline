#!/bin/bash
# Bioinformatics Project Setup Script
# This script sets up the complete environment for bioinformatics analysis

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Script metadata
SCRIPT_NAME="setup.sh"
SCRIPT_VERSION="1.0.0"
SCRIPT_AUTHOR="bioinfo-project"
SCRIPT_DATE="2024-12-16"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

# Banner
print_banner() {
    echo "================================================"
    echo "    Bioinformatics Project Setup Script"
    echo "    Version: $SCRIPT_VERSION"
    echo "    Author: $SCRIPT_AUTHOR"
    echo "    Date: $SCRIPT_DATE"
    echo "================================================"
    echo
}

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Bioinformatics project setup script

OPTIONS:
    -h, --help          Show this help message
    -v, --verbose       Enable verbose output
    -f, --force         Force overwrite existing installations
    --skip-conda        Skip conda environment setup
    --skip-containers   Skip container builds
    --skip-r            Skip R environment setup
    --skip-data         Skip test data download
    --conda-only        Only setup conda environments
    --containers-only   Only build containers
    --check-deps        Only check dependencies
    --clean             Clean existing setup before installing

EXAMPLES:
    $0                  # Full setup
    $0 --conda-only     # Only setup conda environments
    $0 --check-deps     # Check if dependencies are installed
    $0 --clean          # Clean and reinstall everything

ENVIRONMENT VARIABLES:
    PROJ_ROOT          Project root directory (default: current directory)
    CONDA_ENV_PREFIX   Conda environment prefix (default: conda/envs)
    SKIP_CONFIRMATION  Skip confirmation prompts (set to 'yes')

EOF
}

# Default values
VERBOSE=false
FORCE=false
SKIP_CONDA=false
SKIP_CONTAINERS=false
SKIP_R=false
SKIP_DATA=false
CONDA_ONLY=false
CONTAINERS_ONLY=false
CHECK_DEPS_ONLY=false
CLEAN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        --skip-conda)
            SKIP_CONDA=true
            shift
            ;;
        --skip-containers)
            SKIP_CONTAINERS=true
            shift
            ;;
        --skip-r)
            SKIP_R=true
            shift
            ;;
        --skip-data)
            SKIP_DATA=true
            shift
            ;;
        --conda-only)
            CONDA_ONLY=true
            shift
            ;;
        --containers-only)
            CONTAINERS_ONLY=true
            shift
            ;;
        --check-deps)
            CHECK_DEPS_ONLY=true
            shift
            ;;
        --clean)
            CLEAN=true
            shift
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
CONDA_ENV_PREFIX=${CONDA_ENV_PREFIX:-"$PROJ_ROOT/envs/conda"}

# Verbose logging
vlog() {
    if [ "$VERBOSE" = true ]; then
        log_info "[VERBOSE] $1"
    fi
}

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check dependencies
check_dependencies() {
    log_info "Checking system dependencies..."
    
    local missing_deps=()
    
    # Essential tools
    local required_tools=("git" "curl" "wget" "tar" "gzip")
    
    for tool in "${required_tools[@]}"; do
        if ! command_exists "$tool"; then
            missing_deps+=("$tool")
        else
            vlog "✓ $tool found"
        fi
    done
    
    # Optional but recommended tools
    local optional_tools=("conda" "mamba" "docker" "singularity" "nextflow" "snakemake")
    
    for tool in "${optional_tools[@]}"; do
        if command_exists "$tool"; then
            vlog "✓ $tool found"
        else
            log_warning "○ $tool not found (optional)"
        fi
    done
    
    # Check conda/mamba
    if ! command_exists "conda" && ! command_exists "mamba"; then
        log_warning "Neither conda nor mamba found. Consider installing Miniconda or Mambaforge."
    fi
    
    # Check container runtimes
    if ! command_exists "docker" && ! command_exists "singularity" && ! command_exists "podman"; then
        log_warning "No container runtime found. Consider installing Docker, Singularity, or Podman."
    fi
    
    if [ ${#missing_deps[@]} -gt 0 ]; then
        log_error "Missing required dependencies: ${missing_deps[*]}"
        log_error "Please install missing dependencies and try again."
        return 1
    fi
    
    log_success "All required dependencies found!"
    return 0
}

# Setup conda environments
setup_conda_environments() {
    log_info "Setting up conda environments..."
    
    if ! command_exists "conda" && ! command_exists "mamba"; then
        log_error "Conda/mamba not found. Please install Miniconda or Mambaforge first."
        return 1
    fi
    
    # Use mamba if available (faster)
    local conda_cmd="conda"
    if command_exists "mamba"; then
        conda_cmd="mamba"
        log_info "Using mamba for faster environment creation"
    fi
    
    # Environment YAML files
    local env_files=(
        "$PROJ_ROOT/envs/conda/base.yml"
        "$PROJ_ROOT/envs/conda/rnaseq.yml"
        "$PROJ_ROOT/envs/conda/genomics.yml"
        "$PROJ_ROOT/envs/conda/ml.yml"
    )
    
    for env_file in "${env_files[@]}"; do
        if [ ! -f "$env_file" ]; then
            log_error "Environment file not found: $env_file"
            continue
        fi
        
        local env_name
        env_name=$(grep "^name:" "$env_file" | cut -d':' -f2 | tr -d ' ')
        
        log_info "Creating environment: $env_name"
        
        # Check if environment already exists
        if conda env list | grep -q "^$env_name "; then
            if [ "$FORCE" = true ]; then
                log_warning "Removing existing environment: $env_name"
                conda env remove -n "$env_name" -y
            else
                log_warning "Environment $env_name already exists. Use --force to recreate."
                continue
            fi
        fi
        
        vlog "Creating environment from: $env_file"
        if $conda_cmd env create -f "$env_file"; then
            log_success "Environment $env_name created successfully"
        else
            log_error "Failed to create environment: $env_name"
        fi
    done
}

# Setup R environment
setup_r_environment() {
    log_info "Setting up R environment..."
    
    if ! command_exists "R"; then
        log_warning "R not found. Installing via conda base environment..."
        conda activate bioinfo-base 2>/dev/null || true
        conda install -c conda-forge r-base -y
    fi
    
    # Check if renv is available
    if R --slave -e "if (!require('renv', quietly=TRUE)) quit(status=1)" 2>/dev/null; then
        log_info "renv package found"
    else
        log_info "Installing renv package..."
        R --slave -e "install.packages('renv', repos='https://cran.rstudio.com/')"
    fi
    
    # Initialize renv in the R environment directory
    local renv_dir="$PROJ_ROOT/envs/renv"
    if [ -d "$renv_dir" ]; then
        log_info "Initializing renv in $renv_dir"
        cd "$renv_dir"
        
        # Only initialize if not already done
        if [ ! -d "renv" ]; then
            R --slave -e "renv::init(bare=TRUE)"
            log_success "renv initialized"
        else
            log_info "renv already initialized"
        fi
        
        cd "$PROJ_ROOT"
    fi
}

# Build containers
build_containers() {
    log_info "Building container images..."
    
    # Docker containers
    if command_exists "docker"; then
        local dockerfile="$PROJ_ROOT/containers/docker/Dockerfile.full"
        if [ -f "$dockerfile" ]; then
            log_info "Building Docker image..."
            cd "$PROJ_ROOT/containers/docker"
            
            if docker build -t bioinfo:latest -f Dockerfile.full .; then
                log_success "Docker image built successfully"
            else
                log_error "Failed to build Docker image"
            fi
            
            cd "$PROJ_ROOT"
        fi
    else
        log_warning "Docker not found. Skipping Docker image build."
    fi
    
    # Singularity containers
    if command_exists "singularity"; then
        local singularity_def="$PROJ_ROOT/containers/singularity/full.def"
        if [ -f "$singularity_def" ]; then
            log_info "Building Singularity image..."
            cd "$PROJ_ROOT/containers/singularity"
            
            if singularity build bioinfo.sif full.def; then
                log_success "Singularity image built successfully"
            else
                log_error "Failed to build Singularity image"
            fi
            
            cd "$PROJ_ROOT"
        fi
    else
        log_warning "Singularity not found. Skipping Singularity image build."
    fi
}

# Download test data
download_test_data() {
    log_info "Setting up test data..."
    
    local data_dir="$PROJ_ROOT/data"
    
    # Create data directories
    mkdir -p "$data_dir/{raw,processed,external,reference,interim}"
    
    # Create .gitkeep files
    touch "$data_dir/raw/.gitkeep"
    touch "$data_dir/processed/.gitkeep"
    touch "$data_dir/external/.gitkeep"
    
    # Create sample test files (small FASTQ files for testing)
    local test_dir="$data_dir/test"
    mkdir -p "$test_dir"
    
    # Create minimal test FASTQ files
    cat > "$test_dir/sample1_R1.fastq" << 'EOF'
@test_read_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@test_read_2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    
    cat > "$test_dir/sample1_R2.fastq" << 'EOF'
@test_read_1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@test_read_2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    
    # Compress test files
    gzip "$test_dir/sample1_R1.fastq"
    gzip "$test_dir/sample1_R2.fastq"
    
    # Create sample sheet
    cat > "$data_dir/samples.tsv" << 'EOF'
sample_id	r1_path	r2_path	condition
sample1	data/test/sample1_R1.fastq.gz	data/test/sample1_R2.fastq.gz	control
EOF
    
    log_success "Test data created in $test_dir"
}

# Cleanup function
cleanup_existing() {
    log_info "Cleaning existing setup..."
    
    if [ "$SKIP_CONFIRMATION" != "yes" ]; then
        read -p "This will remove existing conda environments and containers. Continue? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "Cleanup cancelled"
            return 0
        fi
    fi
    
    # Remove conda environments
    local envs=("bioinfo-base" "bioinfo-rnaseq" "bioinfo-genomics" "bioinfo-ml")
    for env in "${envs[@]}"; do
        if conda env list | grep -q "^$env "; then
            log_info "Removing conda environment: $env"
            conda env remove -n "$env" -y
        fi
    done
    
    # Remove containers
    if command_exists "docker"; then
        if docker images | grep -q "bioinfo"; then
            log_info "Removing Docker images..."
            docker rmi bioinfo:latest 2>/dev/null || true
        fi
    fi
    
    if [ -f "$PROJ_ROOT/containers/singularity/bioinfo.sif" ]; then
        log_info "Removing Singularity image..."
        rm -f "$PROJ_ROOT/containers/singularity/bioinfo.sif"
    fi
    
    log_success "Cleanup completed"
}

# Main setup function
main_setup() {
    log_info "Starting bioinformatics project setup..."
    log_info "Project root: $PROJ_ROOT"
    
    # Change to project root
    cd "$PROJ_ROOT"
    
    # Check dependencies first
    if ! check_dependencies; then
        return 1
    fi
    
    if [ "$CHECK_DEPS_ONLY" = true ]; then
        log_success "Dependency check completed"
        return 0
    fi
    
    # Clean if requested
    if [ "$CLEAN" = true ]; then
        cleanup_existing
    fi
    
    # Setup based on options
    if [ "$CONTAINERS_ONLY" = true ]; then
        build_containers
    elif [ "$CONDA_ONLY" = true ]; then
        setup_conda_environments
    else
        # Full setup
        if [ "$SKIP_CONDA" = false ]; then
            setup_conda_environments
        fi
        
        if [ "$SKIP_R" = false ]; then
            setup_r_environment
        fi
        
        if [ "$SKIP_CONTAINERS" = false ]; then
            build_containers
        fi
        
        if [ "$SKIP_DATA" = false ]; then
            download_test_data
        fi
    fi
    
    log_success "Setup completed successfully!"
    
    # Print next steps
    echo
    log_info "Next steps:"
    echo "  1. Activate a conda environment: conda activate bioinfo-base"
    echo "  2. Run a test workflow: bash bin/run_analysis.sh"
    echo "  3. Check the documentation in docs/"
    echo "  4. Explore the example data in data/test/"
    echo
}

# Error handling
trap 'log_error "Script failed at line $LINENO"' ERR

# Main execution
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    print_banner
    main_setup
fi

