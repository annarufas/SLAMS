#!/usr/bin/env bash

# ======================================================================================
# SLAMS – SLURM Worker for Balanced Runs
#
# Author: A. Rufas
# Created: 16 Feb 2026    
# Last updated: 16 Sep 2026 
#
# PURPOSE
# -------
#	Executes one SLAMS run directory. This worker supports two modes:
#
#   1. MAPPED ARRAY MODE - normal production batches
#      run_model_balanced.sh submits a job array and supplies MAPPING_FILE.
#      SLURM_ARRAY_TASK_ID is looked up in that file to obtain the actual
#      run number. This supports non-contiguous run IDs.
#
#   2. DIRECT RUN MODE - targeted re-runs
#      check_output_status.sh submits one job and supplies RUN_ID directly.
#      No mapping file is needed in this mode.
#
# NOTES
# -----
#   - This is a low-level execution worker invoked through SLURM.
#   - Each job executes exactly one run_i directory.
#   - Exactly one of MAPPING_FILE or RUN_ID must be supplied.
#
# ======================================================================================

set -euo pipefail

# Determine absolute path to this script.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR/../../.."
else
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    PROJECT_ROOT="$( cd "$SCRIPT_DIR/../../.." && pwd )"
fi

CORE_DIR="$PROJECT_ROOT/code/bash/core"

# Source. 
source "$CORE_DIR/slams_logging.sh"
source "$CORE_DIR/slams_paths_config.sh"
source "$CORE_DIR/slams_workflow_config.sh"
source "$CORE_DIR/slams_environment_init.sh"
slams_paths_config 
slams_workflow_config
slams_environment_init
	
# Validate variables shared by both execution modes.
: "${CONFIG_NAME:?Missing CONFIG_NAME}"
: "${RUNS_DIR:?Missing RUNS_DIR}"
: "${INTEL_OMP_LIB:?Missing INTEL_OMP_LIB}"

# Resolve the actual run ID according to the submission mode.
if [[ -n "${MAPPING_FILE:-}" && -n "${RUN_ID:-}" ]]; then
    die "Both MAPPING_FILE and RUN_ID were supplied; provide only one"

elif [[ -n "${MAPPING_FILE:-}" ]]; then
    # Mapped array mode:
    # SLURM_ARRAY_TASK_ID is a lookup key. It may not equal the actual
    # run number because run IDs can be non-contiguous.
    : "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID in mapped array mode}"

    RUN_ID=$(awk -v slurm_id="$SLURM_ARRAY_TASK_ID" '$1 == slurm_id {print $2}' "$MAPPING_FILE")

    if [[ -z "$RUN_ID" ]]; then
        die "Could not find RUN_ID for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID in $MAPPING_FILE"
    fi
    log "Mapped array mode: array task $SLURM_ARRAY_TASK_ID corresponds to run_$RUN_ID"

elif [[ -n "${RUN_ID:-}" ]]; then
    # Direct run mode:
    # RUN_ID was supplied directly by check_output_status.sh.
    log "Direct run mode: executing run_$RUN_ID"

else
    die "Either MAPPING_FILE or RUN_ID must be supplied"
fi

# Define the directory for this specific run.
thisRunDir="$RUNS_DIR"/run_"$RUN_ID"
if [[ ! -d "$thisRunDir" ]]; then
    die "Run directory $thisRunDir does not exist"
fi

# Change to the directory.
cd "$thisRunDir" 
log "Running job ID $RUN_ID in $(pwd)"

# Clean old outputs.
rm -f ./*.out ./log* ./out*

# Execute SLAMS model.
# SLURM captures stdout/stderr using the --output option from the submitting job.
srun --export=ALL,LD_LIBRARY_PATH="${INTEL_OMP_LIB}:${LD_LIBRARY_PATH:-}" \
     -n 1 ./SLAMSexecutable