#!/usr/bin/env bash

# ======================================================================================
# SLAMS 2.0 – SLURM Job-Array Worker (Balanced Runs)
#
# Author: A. Rufas | 16 Feb 2026
#
# PURPOSE
# -------
#   Executes a single SLAMS run directory as part of a SLURM job array.
#     - Resolve the correct run_* directory using a mapping file
#     - Validate required environment variables
#     - Cleans previous model outputs
#     - Executes SLAMSexecutable using srun
#
# ARCHITECTURAL ROLE
# ------------------
#   This is a low-level execution worker.
#   It is NOT part of the canonical pipeline.
#   It is invoked exclusively by run_model_balanced.sh via SLURM.
#
#   Each SLURM array task executes exactly one run_i directory.
#
# ======================================================================================

set -euo pipefail

# Determine absolute path to this script
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR/../../.."
else
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    PROJECT_ROOT="$( cd "$SCRIPT_DIR/../../.." && pwd )"
fi

CORE_DIR="$PROJECT_ROOT/code/bash/core"

# Source 
source "$CORE_DIR/slams_logging.sh"
source "$CORE_DIR/slams_paths_config.sh"
source "$CORE_DIR/slams_workflow_config.sh"
source "$CORE_DIR/slams_environment_init.sh"
slams_paths_config 
slams_workflow_config
slams_environment_init
	
# Validate required env vars
: "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
: "${MAPPING_FILE:?Missing MAPPING_FILE}"
: "${RUNS_DIR:?Missing RUNS_DIR}"
: "${INTEL_OMP_LIB:?Missing INTEL_OMP_LIB}"

# Resolve actual run ID from mapping file.
# SLURM_ARRAY_TASK_ID indexes the array (1..N), but actual run IDs may be non-contiguous.
RUN_ID=$(awk -v slurm_id="$SLURM_ARRAY_TASK_ID" '$1 == slurm_id {print $2}' "$MAPPING_FILE")
if [[ -z "$RUN_ID" ]]; then
    die "Could not find RUN_ID for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID in $MAPPING_FILE"
fi

# Define the directory for this specific run
thisRunDir="$RUNS_DIR"/run_"$RUN_ID"
if [[ ! -d "$thisRunDir" ]]; then
    die "Run directory $thisRunDir does not exist"
fi

# Change to the directory
cd "$thisRunDir" 
log "Running job ID $RUN_ID in $(pwd)"

# Clean old outputs
make cleanoutput || true

# Execute SLAMS model.
# SLURM captures stdout/stderr via job array logging.
srun --export=ALL,LD_LIBRARY_PATH="${INTEL_OMP_LIB}:${LD_LIBRARY_PATH}" \
     -n 1 ./SLAMSexecutable