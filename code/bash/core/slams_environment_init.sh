#!/usr/bin/env bash

# ====================================================================================== 
# SLAMS – Environment Initialisation
#
# Author: A. Rufas
# Created: 16 Feb 2026    
# Last updated: 14 Sep 2026 
#                                                                                                                                                                           
# PURPOSE
# -------
#   Initialises the runtime environment for SLAMS by:
#     - Validating CONFIG_NAME
#     - Loading required modules (MATLAB, Intel, MKL) when running under SLURM
#     - Configuring library paths (e.g. MKL, OpenMP)
#
# NOTES
# -----
# 	- This file contains internal functions used by SLAMS workflows.
# 	- It is sourced by higher-level pipeline scripts and is not intended to be executed or 
# 	modified directly by users.
# 	- HPC module names and versions are site-specific and may need to be adapted for 
#   different systems.
#   
# ======================================================================================                

slams_environment_init () 
{	
    if [[ -n "${SLURM_JOB_ID:-}" || "${CHOICE_RUN_IN_SLURM:-false}" == "true" ]]; then
        log "Loading modules for SLURM execution"

        # --- MATLAB --- 
        module load MATLAB/R2022a || {
            die "Failed to load MATLAB"
        }
        
        # --- Intel compiler (first match wins) ---
        INTEL_MODULE=""
        for m in intel-compilers/2021.4.0 intel/2022a intel/2022.1 intel; do
            if module load "$m"; then
                INTEL_MODULE="$m"
                log "Loaded $m"
                break
            fi
        done
        if [[ -z "$INTEL_MODULE" ]]; then
            die "No Intel compiler module loaded"
        fi

        # --- Intel MKL (first match wins) ---
        for m in imkl/2021.4.0 imkl/2022.1.0 imkl/2022.2.1 imkl; do
            module load "$m" && { log "Loaded $m"; break; }
        done

        # --- Canonical MKL configuration ---
        MKL_VERSION="2021.4.0"
        export MKL_VERSION
        if [[ -z "${EBROOTIMKL:-}" ]]; then
            die "EBROOTIMKL not set after loading imkl module"
        fi
        MKL_LIBDIR="${EBROOTIMKL}/mkl/${MKL_VERSION}/lib/intel64"
        [[ -d "$MKL_LIBDIR" ]] || {
            die "MKL_LIBDIR not found: $MKL_LIBDIR"
        }

        export MKL_LIBDIR
        export LD_LIBRARY_PATH="${MKL_LIBDIR}:${LD_LIBRARY_PATH:-}"

        # --- Intel OpenMP runtime (libiomp5.so) ---
        INTEL_OMP_LIB=$(find "${EBROOTINTELMINCOMPILERS}/compiler" \
    		-path "*/linux/compiler/lib/intel64_lin" -type d | head -n1)

		if [[ -z "${INTEL_OMP_LIB:-}" || ! -f "${INTEL_OMP_LIB}/libiomp5.so" ]]; then
    		die "libiomp5.so not found in Intel compiler path"
		fi

		export INTEL_OMP_LIB
		export LD_LIBRARY_PATH="${INTEL_OMP_LIB}:${LD_LIBRARY_PATH:-}"
    
    else
    	mkdir -p "$LOGS_DIR"
    	exec > "$LOGS_DIR/log_${CONFIG_NAME}.log" 2>&1
    fi
	
	log "Environment initialised successfully"
}