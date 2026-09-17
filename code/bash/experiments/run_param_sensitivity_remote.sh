#!/usr/bin/env bash
#SBATCH --job-name=sensit 
#SBATCH --time=00:10:00                   
#SBATCH --partition=devel
#SBATCH --output=../../../logs/log_sensit_%A.log

set -euo pipefail # abort the script if errors with commands, variables and pipelines occur
IFS=$'\n\t'       # controls bash word splitting

# ======================================================================================
# SLAMS – Sensitivity Analysis Job
#
# Author: A. Rufas
# Created: 16 Feb 2026
# Last updated: 14 Sep 2026
#
# PURPOSE
# -------
# 	This script performs a parameter sensitivity analysis by:
#   	1. Modifying a base namelist using values specified in CHANGES_FILE
#   	2. Creating corresponding latitude/longitude configuration files
#   	3. Submitting SLURM jobs for model input creation, model execution, and post-processing
#      
# IMPORTANT: EXECUTION MODEL
# --------------------------
# 	This script is NOT intended to be run end-to-end in a single submission. Several steps
# 	must be executed sequentially, with downstream jobs submitted only after upstream jobs 
# 	have completed successfully.
#
# 	The execution order is controlled by commenting/uncommenting function calls inside the
# 	main() function.
#
# REQUIRED PREPARATION
# --------------------
# 	1. In slams_workflow_config.sh:
#    	- Ensure CHOICE_RECOMPILE=false so the model is not recompiled for every 
#	s	ensitivity run.
#      
# 	2. Verify that CHANGES_FILE exists in config/ and contains:
#      	parameter,value,bound
#    	as a CSV with a header row.
#
# 	3. In the template namelist:
#    	- Ensure CHOICE_GRID_DOMAIN=2 so the model is run on a local basis.
#
# STEP-BY-STEP WORKFLOW
# --------------------
# 	Step 1: Generate namelists and coordinate files
#   	- Uncomment ONLY the following in main():
#       	create_new_namelist
#       	create_corresponding_coords_file
#   	- Run this script and verify output files.
#
# 	Step 2: Create model input (must finish before Step 3)
#   	- Uncomment:
#       	submit_script_create_model_input
#   	- Comment out all other submit_* functions
#   	- Submit the script and wait until ALL jobs finish.
#
# 	Step 3: Run the balanced model
#   	- Uncomment:
#       	submit_script_run_model_balanced
#   	- Submit the script and wait until ALL jobs finish.
#
# 	Step 4: Check output status
#   	- Uncomment:
#       	submit_check_output_status
#
# 	Step 5: Read and post-process model output
#   	- Uncomment:
#       	submit_read_model_output
# 
# 	Step 6: Download model output to local machine
#   	- Use the script download_runsoutput.sh locally
#
# SLURM NOTES
# -----------
# 	- No memory directive (#SBATCH --mem) is specified intentionally.
# 	- This prevents memory settings from being inherited by downstream jobs.
# 	- This script itself performs only lightweight Bash operations.
#
# USAGE
# -----
# 	On a local machine (for file generation only):
#   	./run_param_sensitivity_test.sh
#
# 	On a SLURM system:
#  		sbatch run_param_sensitivity_test.sh
#
# ======================================================================================

# --------------------------------------------------
# Configuration
# --------------------------------------------------

# Determine absolute path to this script
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR/../../.."
else
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    PROJECT_ROOT="$( cd "$SCRIPT_DIR/../../.." && pwd )"
fi

CORE_DIR="$PROJECT_ROOT/code/bash/core"
PIPELINE_DIR="$PROJECT_ROOT/code/bash/pipelines"
MAINTENANCE_DIR="$PROJECT_ROOT/code/bash/maintenance"
CONFIG_DIR="$PROJECT_ROOT/config"
BASE_CONFIG_NAME="LOCALTS6"
CHANGES_FILE="$CONFIG_DIR/parameter_sensitivity_analysis.csv"
export CONFIG_DIR BASE_CONFIG_NAME CHANGES_FILE
source "$CORE_DIR/slams_logging.sh"

sed -i -e '$a\' "$CHANGES_FILE" # ensure last line has newline

# --------------------------------------------------
# Establish run order (uncomment as needed)
# --------------------------------------------------

main () 
{
	#create_new_namelist
	#create_corresponding_coords_file
	#submit_script_create_model_input
	#submit_script_run_model_balanced
	#submit_check_output_status
	submit_read_model_output
}

# --------------------------------------------------
# Functions used in this script
# --------------------------------------------------

create_new_namelist () 
{
	local base_namelist="$CONFIG_DIR/namelist_${BASE_CONFIG_NAME}.txt"

	if [[ ! -f "$CHANGES_FILE" ]]; then
    	die "Changes file not found at $CHANGES_FILE" >&2
	fi
	
	# Read each line from the changes file, skipping the header
	tail -n +2 "$CHANGES_FILE" | while IFS=',' read -r param value bound
	do
		# Strip surrounding quotes and carriage returns (e.g. from Windows line endings)
    	param=$(echo "$param" | tr -d '"' | tr -d '\r')
    	value=$(echo "$value" | tr -d '"' | tr -d '\r')
    	bound=$(echo "$bound" | tr -d '"' | tr -d '\r')
    
		local output_file="$CONFIG_DIR/namelist_${BASE_CONFIG_NAME}_${param}_${bound}.txt"
        local escaped_param
        escaped_param=$(printf '%s\n' "$param" | sed -e 's/[]\/$*.^[]/\\&/g')

        sed -E "s/^(${escaped_param}[[:space:]]*=[[:space:]]*)[^[:space:]]+/\1${value}/" "$base_namelist" > "$output_file"
        log "Created: $output_file"
	done
	log "New namelist created"
}

create_corresponding_coords_file ()
{
    local base_lat_coords_file="$CONFIG_DIR/config_latitudes_${BASE_CONFIG_NAME}.txt"
    local base_lon_coords_file="$CONFIG_DIR/config_longitudes_${BASE_CONFIG_NAME}.txt"

    # Read each line from the changes file, skipping the header
	tail -n +2 "$CHANGES_FILE" | while IFS=',' read -r param value bound
    do
    	# Strip surrounding quotes and carriage returns (e.g. from Windows line endings)
    	param=$(echo "$param" | tr -d '"' | tr -d '\r')
    	value=$(echo "$value" | tr -d '"' | tr -d '\r')
    	bound=$(echo "$bound" | tr -d '"' | tr -d '\r')
    
        local output_lat_coords_file="$CONFIG_DIR/config_latitudes_${BASE_CONFIG_NAME}_${param}_${bound}.txt"
        cp "$base_lat_coords_file" "$output_lat_coords_file"

        local output_lon_coords_file="$CONFIG_DIR/config_longitudes_${BASE_CONFIG_NAME}_${param}_${bound}.txt"
        cp "$base_lon_coords_file" "$output_lon_coords_file"
    done
    log "New coords file created"
}

submit_script_create_model_input ()
{
    # Read each line from the changes file, skipping the header
	tail -n +2 "$CHANGES_FILE" | while IFS=',' read -r param value bound
    do
    	# Strip surrounding quotes and carriage returns (e.g. from Windows line endings)
    	param=$(echo "$param" | tr -d '"' | tr -d '\r')
    	value=$(echo "$value" | tr -d '"' | tr -d '\r')
    	bound=$(echo "$bound" | tr -d '"' | tr -d '\r')
    	
        local new_config_name="${BASE_CONFIG_NAME}_${param}_${bound}"
        jobid=$(sbatch --export=ALL,CONFIG_NAME="$new_config_name" "$PIPELINE_DIR/create_model_input.sh" | awk '{print $NF}')
        log "Submitted: $new_config_name (Job ID: $jobid)"
    done
    log "Submitting script that creates model input"
}

submit_script_run_model_balanced ()
{
    # Read each line from the changes file, skipping the header
	tail -n +2 "$CHANGES_FILE" | while IFS=',' read -r param value bound
    do
    	# Strip surrounding quotes and carriage returns (e.g. from Windows line endings)
    	param=$(echo "$param" | tr -d '"' | tr -d '\r')
    	value=$(echo "$value" | tr -d '"' | tr -d '\r')
    	bound=$(echo "$bound" | tr -d '"' | tr -d '\r')
    	
        local new_config_name="${BASE_CONFIG_NAME}_${param}_${bound}"
        jobid=$(sbatch --export=ALL,CONFIG_NAME="$new_config_name" "$PIPELINE_DIR/run_model_balanced.sh" | awk '{print $NF}')
        log "Submitted: $new_config_name (Job ID: $jobid)"
    done
    log "Runs submitted"
}

submit_check_output_status ()
{
    # Read each line from the changes file, skipping the header
	tail -n +2 "$CHANGES_FILE" | while IFS=',' read -r param value bound
    do
    	# Strip surrounding quotes and carriage returns (e.g. from Windows line endings)
    	param=$(echo "$param" | tr -d '"' | tr -d '\r')
    	value=$(echo "$value" | tr -d '"' | tr -d '\r')
    	bound=$(echo "$bound" | tr -d '"' | tr -d '\r')
    	
        local new_config_name="${BASE_CONFIG_NAME}_${param}_${bound}"
        jobid=$(sbatch --export=ALL,CONFIG_NAME="$new_config_name" "$MAINTENANCE_DIR/check_output_status.sh" | awk '{print $NF}')
        log "Submitted: $new_config_name (Job ID: $jobid)"
    done
}

submit_read_model_output ()
{
    # Read each line from the changes file, skipping the header
	tail -n +2 "$CHANGES_FILE" | while IFS=',' read -r param value bound
    do
    	# Strip surrounding quotes and carriage returns (e.g. from Windows line endings)
    	param=$(echo "$param" | tr -d '"' | tr -d '\r')
    	value=$(echo "$value" | tr -d '"' | tr -d '\r')
    	bound=$(echo "$bound" | tr -d '"' | tr -d '\r')
    	
        local new_config_name="${BASE_CONFIG_NAME}_${param}_${bound}"
        jobid=$(sbatch --export=ALL,CONFIG_NAME="$new_config_name" "$PIPELINE_DIR/read_model_output.sh" | awk '{print $NF}')
        log "Submitted: $new_config_name (Job ID: $jobid)"
    done
}

# Execute this script
main