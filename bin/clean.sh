#!/bin/bash
# Bioinformatics Project Cleanup Script
# This script cleans up temporary files and outputs

set -euo pipefail

# Script metadata
SCRIPT_NAME="clean.sh"
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
CLEAN_RESULTS=false
CLEAN_LOGS=false
CLEAN_CACHE=false
CLEAN_TEMP=true
CLEAN_CONTAINERS=false
CLEAN_ENVS=false
DRY_RUN=false
FORCE=false
VERBOSE=false

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Clean up bioinformatics project files and directories

OPTIONS:
    -h, --help              Show this help message
    -r, --results           Clean results directory
    -l, --logs              Clean log files
    -c, --cache             Clean cache directories
    -t, --temp              Clean temporary files [default: enabled]
    --containers            Clean container images
    --envs                  Clean conda environments
    -a, --all               Clean everything (use with caution)
    -n, --dry-run          Show what would be cleaned without doing it
    -f, --force            Force cleanup without confirmation
    -v, --verbose          Enable verbose output

EXAMPLES:
    $0                      # Clean temporary files only
    $0 -t -l                # Clean temp files and logs
    $0 --all --dry-run      # Show what would be cleaned
    $0 --results --force    # Force clean results directory

WARNING:
    Using --all or --results will permanently delete analysis outputs!
    Always use --dry-run first to see what will be deleted.

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -r|--results)
            CLEAN_RESULTS=true
            shift
            ;;
        -l|--logs)
            CLEAN_LOGS=true
            shift
            ;;
        -c|--cache)
            CLEAN_CACHE=true
            shift
            ;;
        -t|--temp)
            CLEAN_TEMP=true
            shift
            ;;
        --containers)
            CLEAN_CONTAINERS=true
            shift
            ;;
        --envs)
            CLEAN_ENVS=true
            shift
            ;;
        -a|--all)
            CLEAN_RESULTS=true
            CLEAN_LOGS=true
            CLEAN_CACHE=true
            CLEAN_TEMP=true
            CLEAN_CONTAINERS=true
            shift
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
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
cd "$PROJ_ROOT"

# Verbose logging
vlog() {
    if [[ "$VERBOSE" = true ]]; then
        log_info "[VERBOSE] $1"
    fi
}

# Dry run logging
dry_run_log() {
    if [[ "$DRY_RUN" = true ]]; then
        log_info "[DRY RUN] Would $1"
    else
        log_info "$1"
    fi
}

# Check if directory exists and is not empty
dir_exists_and_not_empty() {
    [[ -d "$1" ]] && [[ $(find "$1" -mindepth 1 -print -quit 2>/dev/null) ]]
}

# Get directory size
get_dir_size() {
    if [[ -d "$1" ]]; then
        du -sh "$1" 2>/dev/null | cut -f1 || echo "unknown"
    else
        echo "0"
    fi
}

# Clean temporary files
clean_temp_files() {
    log_info "Cleaning temporary files..."
    
    local temp_patterns=(
        "tmp/*"
        "temp/*"
        "cache/tmp/*"
        ".nextflow/*"
        ".nextflow.log*"
        ".snakemake/*"
        "work/*"
        "*.tmp"
        "*~"
        "*.swp"
        "*.swo"
        "*.bak"
        "*.orig"
        "*_temp"
        "*_tmp"
    )
    
    local total_cleaned=0
    
    for pattern in "${temp_patterns[@]}"; do
        if ls $pattern 1> /dev/null 2>&1; then
            local size
            size=$(du -sch $pattern 2>/dev/null | tail -1 | cut -f1 || echo "unknown")
            
            dry_run_log "remove temp files matching: $pattern ($size)"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -rf $pattern
                ((total_cleaned++))
            fi
        fi
    done
    
    # Clean workflow work directories
    local work_dirs=("work" "tmp" "temp")
    for dir in "${work_dirs[@]}"; do
        if dir_exists_and_not_empty "$dir"; then
            local size
            size=$(get_dir_size "$dir")
            
            dry_run_log "remove work directory: $dir ($size)"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -rf "$dir"
                ((total_cleaned++))
            fi
        fi
    done
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Cleaned $total_cleaned temporary file groups"
    fi
}

# Clean log files
clean_logs() {
    log_info "Cleaning log files..."
    
    local log_dirs=("logs" ".nextflow/log")
    local log_patterns=("*.log" "*.out" "*.err")
    
    for dir in "${log_dirs[@]}"; do
        if dir_exists_and_not_empty "$dir"; then
            local size
            size=$(get_dir_size "$dir")
            
            dry_run_log "remove log directory: $dir ($size)"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -rf "$dir"/*
            fi
        fi
    done
    
    # Clean individual log files
    for pattern in "${log_patterns[@]}"; do
        if ls $pattern 1> /dev/null 2>&1; then
            dry_run_log "remove log files: $pattern"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -f $pattern
            fi
        fi
    done
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Log files cleaned"
    fi
}

# Clean cache directories
clean_cache() {
    log_info "Cleaning cache directories..."
    
    local cache_dirs=(
        "cache"
        ".cache"
        ".conda/pkgs"
        ".singularity"
        ".nextflow/cache"
        "__pycache__"
    )
    
    for dir in "${cache_dirs[@]}"; do
        if dir_exists_and_not_empty "$dir"; then
            local size
            size=$(get_dir_size "$dir")
            
            dry_run_log "remove cache directory: $dir ($size)"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -rf "$dir"
            fi
        fi
    done
    
    # Clean Python cache files
    find . -name "__pycache__" -type d 2>/dev/null | while read -r dir; do
        dry_run_log "remove Python cache: $dir"
        
        if [[ "$DRY_RUN" = false ]]; then
            rm -rf "$dir"
        fi
    done
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Cache directories cleaned"
    fi
}

# Clean results
clean_results() {
    log_warning "Cleaning results directories - this will delete analysis outputs!"
    
    local results_dirs=("results" "output" "outputs")
    
    if [[ "$FORCE" = false ]] && [[ "$DRY_RUN" = false ]]; then
        read -p "Are you sure you want to delete results? This cannot be undone! (type 'yes'): " -r
        if [[ $REPLY != "yes" ]]; then
            log_info "Results cleanup cancelled"
            return 0
        fi
    fi
    
    for dir in "${results_dirs[@]}"; do
        if dir_exists_and_not_empty "$dir"; then
            local size
            size=$(get_dir_size "$dir")
            
            dry_run_log "remove results directory: $dir ($size)"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -rf "$dir"
            fi
        fi
    done
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Results directories cleaned"
    fi
}

# Clean containers
clean_containers() {
    log_info "Cleaning container images..."
    
    # Docker images
    if command -v docker >/dev/null 2>&1; then
        local docker_images
        docker_images=$(docker images --filter "reference=bioinfo*" -q 2>/dev/null || true)
        
        if [[ -n "$docker_images" ]]; then
            dry_run_log "remove Docker images: bioinfo*"
            
            if [[ "$DRY_RUN" = false ]]; then
                echo "$docker_images" | xargs docker rmi -f 2>/dev/null || true
            fi
        fi
    fi
    
    # Singularity images
    local sif_files
    sif_files=$(find . -name "*.sif" 2>/dev/null || true)
    
    if [[ -n "$sif_files" ]]; then
        while IFS= read -r sif_file; do
            local size
            size=$(du -sh "$sif_file" 2>/dev/null | cut -f1 || echo "unknown")
            
            dry_run_log "remove Singularity image: $sif_file ($size)"
            
            if [[ "$DRY_RUN" = false ]]; then
                rm -f "$sif_file"
            fi
        done <<< "$sif_files"
    fi
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Container images cleaned"
    fi
}

# Clean conda environments
clean_environments() {
    log_info "Cleaning conda environments..."
    
    if ! command -v conda >/dev/null 2>&1; then
        log_warning "Conda not found, skipping environment cleanup"
        return 0
    fi
    
    local bioinfo_envs
    bioinfo_envs=$(conda env list | grep "bioinfo-" | awk '{print $1}' || true)
    
    if [[ -n "$bioinfo_envs" ]]; then
        while IFS= read -r env_name; do
            if [[ -n "$env_name" ]]; then
                dry_run_log "remove conda environment: $env_name"
                
                if [[ "$DRY_RUN" = false ]]; then
                    conda env remove -n "$env_name" -y 2>/dev/null || true
                fi
            fi
        done <<< "$bioinfo_envs"
    fi
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Conda environments cleaned"
    fi
}

# Show cleanup summary
show_summary() {
    log_info "Cleanup summary:"
    echo "  Temporary files: $([ "$CLEAN_TEMP" = true ] && echo "✓" || echo "○")"
    echo "  Log files: $([ "$CLEAN_LOGS" = true ] && echo "✓" || echo "○")"
    echo "  Cache directories: $([ "$CLEAN_CACHE" = true ] && echo "✓" || echo "○")"
    echo "  Results: $([ "$CLEAN_RESULTS" = true ] && echo "✓" || echo "○")"
    echo "  Containers: $([ "$CLEAN_CONTAINERS" = true ] && echo "✓" || echo "○")"
    echo "  Environments: $([ "$CLEAN_ENVS" = true ] && echo "✓" || echo "○")"
    
    if [[ "$DRY_RUN" = true ]]; then
        echo
        log_info "This was a dry run. Use without --dry-run to actually clean."
    fi
}

# Main execution
main() {
    log_info "Starting cleanup process..."
    log_info "Project root: $PROJ_ROOT"
    
    if [[ "$DRY_RUN" = true ]]; then
        log_info "DRY RUN MODE - No files will be deleted"
    fi
    
    # If no specific options, just clean temp
    if [[ "$CLEAN_RESULTS" = false ]] && [[ "$CLEAN_LOGS" = false ]] && \
       [[ "$CLEAN_CACHE" = false ]] && [[ "$CLEAN_CONTAINERS" = false ]] && \
       [[ "$CLEAN_ENVS" = false ]]; then
        CLEAN_TEMP=true
    fi
    
    # Execute cleanup operations
    if [[ "$CLEAN_TEMP" = true ]]; then
        clean_temp_files
    fi
    
    if [[ "$CLEAN_LOGS" = true ]]; then
        clean_logs
    fi
    
    if [[ "$CLEAN_CACHE" = true ]]; then
        clean_cache
    fi
    
    if [[ "$CLEAN_RESULTS" = true ]]; then
        clean_results
    fi
    
    if [[ "$CLEAN_CONTAINERS" = true ]]; then
        clean_containers
    fi
    
    if [[ "$CLEAN_ENVS" = true ]]; then
        clean_environments
    fi
    
    show_summary
    
    if [[ "$DRY_RUN" = false ]]; then
        log_success "Cleanup completed successfully!"
    fi
}

# Error handling
trap 'log_error "Cleanup failed at line $LINENO"' ERR

# Execute main function
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi

