#!/usr/bin/env bash

#SBATCH --job-name=readout 
#SBATCH --time=14:00:00                   
#SBATCH --partition=long
#SBATCH --qos=earth
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB  # Alternative: total memory instead of per CPU
#SBATCH --output=../../../logs/log_readout_%A.log

set -euo pipefail # abort the script if errors with commands, variables and pipelines occur
IFS=$'\n\t'       # controls bash word splitting

# ======================================================================================
# SLAMS 2.0 – Post-Processing Pipeline Stage
#
# Author: A. Rufas | 16 Feb 2026
#
# PURPOSE
# -------
#   Post-processes SLAMS simulation outputs using MATLAB.
#     - Read model output files from all run_* directories
#     - Aggregates results into an output MATLAB array
#     - Optionally computes derived BCP (Biological Carbon Pump) metrics
#
# USAGE
# -----
#   Local:
#       ./read_model_output.sh <config_name> 
#
#   SLURM:
#       sbatch --export=ALL,CONFIG_NAME=<config_name> read_model_output.sh
# 
# ASSUMPTIONS
# ------------------
#     - All run_* directories exist
#     - Model output files are present
#
# PERFORMANCE NOTES
# -----------------
#   Global configuration:
#       readSlamsOutput:                    ~14 hours
#       calculateBcpMetricsFromSlamsOutput: ~60 hours
#   Local configuration:                    < 5 minutes
#
# ======================================================================================

# --------------------------------------------------
# Logging utilities
# --------------------------------------------------
# log  → standard informational message
# warn → non-fatal issue
# die  → fatal error (exits script)

log()   { printf "[%s] INFO  %s\n"  "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
warn()  { printf "[%s] WARN  %s\n"  "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()   { printf "[%s] ERROR %s\n"  "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

# --------------------------------------------------
# Load core workflow modules
# --------------------------------------------------

# Determine absolute path to this script
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR/../../.."
else
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    PROJECT_ROOT="$( cd "$SCRIPT_DIR/../../.." && pwd )"
fi

CORE_DIR="$PROJECT_ROOT/code/bash/core"

source "$CORE_DIR/slams_logging.sh"
source "$CORE_DIR/slams_paths_config.sh"
source "$CORE_DIR/slams_workflow_config.sh"
source "$CORE_DIR/slams_environment_init.sh"
source "$CORE_DIR/slams_namelist_to_shell.sh"

# --------------------------------------------------
# Validate configuration
# --------------------------------------------------

validate_config_name () 
{
	if [[ -z "${CONFIG_NAME:-}" ]]; then
    	die "CONFIG_NAME is not set"
	fi
	if [[ ! "$CONFIG_NAME" =~ ^[A-Za-z0-9_-]+$ ]]; then
    	die "CONFIG_NAME contains invalid characters"
	fi
}

# --------------------------------------------------
# MATLAB: Read and aggregate SLAMS outputs
# --------------------------------------------------

read_output_with_matlab ()
{
	log "Reading output files, processing the output data and saving them in one big array in $RUNS_DIR..."
	matlab -nodisplay -nodesktop -nosplash -r \
		"addpath(genpath('$MATLABCODE_DIR/')); \
		addpath('$RESOURCES_INTERNAL_DIR/'); \
		addpath(genpath('$RESOURCES_EXTERNAL_DIR/')); \
		addpath('$PROCESSED_DATA_DIR/'); \
		readSlamsOutput('$RUNS_DIR/','$INPUTDATA_DIR/','$FILENAME_SLAMS_OUTPUT',\
					    '$FILENAME_RUN_GRID','$FILENAME_NUM_DEPTHLAYERS',$CHOICE_GRID_DOMAIN,\
						$CHOICE_TEST_MODEL_EQUILIBRIUM,$CHOICE_SAVE_DISK_SPACE_FOR_OUTPUT,\
						$CHOICE_SAVE_SEASONAL_OUTPUT); \
		clc; exit;"
	if [ $? -eq 0 ]; then
		log "Model output files read successfully"
	else
		die "Reading model output files failed"
	fi
}

# --------------------------------------------------
# MATLAB: Compute derived BCP metrics
# --------------------------------------------------

calculate_bcp_metrics_with_matlab ()
{
	export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
	log "Calculating BCP metrics from SLAMS output and saving them in one big array in $RUNS_DIR..."
	matlab -nodisplay -nodesktop -nosplash -r \
		"addpath(genpath('$MATLABCODE_DIR/')); \
		addpath('$RESOURCES_INTERNAL_DIR/'); \
		addpath(genpath('$RESOURCES_EXTERNAL_DIR/')); \
		addpath('$RAW_DATA_DIR/'); \
		calculateBcpMetricsFromSlamsOutput('$RUNS_DIR/','$INPUTDATA_DIR/','$PROCESSED_DATA_DIR/',\
									       '$FILENAME_SLAMS_OUTPUT','$FILENAME_RUN_GRID','$FILENAME_NUM_DEPTHLAYERS',\
									       '$FILENAME_DATA_NPP','$FILENAME_DATA_MASK','$FILENAME_DATA_ZEU',\
										   '$FILENAME_SLAMS_BCPMETRICS',$CHOICE_GRID_DOMAIN); \
    	clc; exit;"
	if [ $? -eq 0 ]; then
		log "BCP metrics calculated successfully"
	else
		die "Calcualting BCP metrics failed"
	fi
}

# --------------------------------------------------
# Pipeline execution order
# --------------------------------------------------
# Defines ordered post-processing steps.
# Additional analysis stages can be enabled as needed.

main () 
{
	log "Starting post-processing stage"
	
	validate_config_name
	
	slams_paths_config
    slams_workflow_config
    
	slams_environment_init
	slams_namelist_to_shell "$CONFIG_NAME" "WorkflowOptions" "ForcingDataFiles" "WorkflowFilenames"

	read_output_with_matlab	
	#calculate_bcp_metrics_with_matlab
	
 	log "Finished"
}

# Execute this script
main