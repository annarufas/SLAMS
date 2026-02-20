#!/usr/bin/env bash

#SBATCH --job-name=runslams 
#SBATCH --time=00:10:00                   
#SBATCH --partition=devel
#SBATCH --output=../../../logs/log_runslams_%A.log

set -euo pipefail # abort the script if errors with commands, variables and pipelines occur
IFS=$'\n\t'       # controls bash word splitting

# ====================================================================================== 
# SLAMS 2.0 – Model Execution Driver (Balanced Runs)
#
# Author: A. Rufas | 16 Feb 2026    
#  
# PURPOSE
# -------
#   Executes a prepared SLAMS experiment across all run_* directories.
#     - Classify runs into short/long categories based on vertical depth (threshold: 200 depth layers)
#     - Split jobs into SLURM arrays for efficient scheduling
#     - Submit job arrays with appropriate memory and time limits
#     - Support limited local execution for small test cases                                                                              
#                                                                                          
# USAGE
# -----
#   Local:
#       ./run_model_balanced.sh <config_name> 
#
#   SLURM:
#       sbatch --export=ALL,CONFIG_NAME=<config_name> run_model_balanced.sh
# 
# ASSUMPTIONS
# -----------
#   - create_model_input.sh has already completed successfully
#   - run_* directories exist in RUNS_DIR
#   - Per-run namelist.input files are present
#   - FILENAME_NUM_DEPTHLAYERS exists in INPUTDATA_DIR
# 
# NOTES
# -----
#   - Partition names, memory, and time limits are defined in manage_runs()
#                                                                                 
# ====================================================================================== 

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
PIPELINE_DIR="$PROJECT_ROOT/code/bash/pipelines"

source "$CORE_DIR/slams_logging.sh"
source "$CORE_DIR/slams_paths_config.sh"
source "$CORE_DIR/slams_workflow_config.sh"
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
# Determine the number of ocean locations
# --------------------------------------------------
# Count number of run_* directories created by forcing generation

count_number_location_directories ()
{
	nLocDirectories=$(find "$INPUTDATA_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
    if [ "$nLocDirectories" -eq 0 ]; then
        die "No location-specific directories found"
    else
		log "Detected $nLocDirectories location-specific directories"
	fi
}

# --------------------------------------------------
# Determine timing of jobs
# --------------------------------------------------

determine_timing_of_jobs_in_slurm ()
{
	if [[ "$CHOICE_RUN_IN_SLURM" == "true" ]]; then

		batchSizeSlurm=1000
		shortJobsCount=0
		longJobsCount=0
		shortBatch=1
		longBatch=1
		
		# Cleanup existing files
		rm -f "$RUNS_DIR"/shortJobs_*.txt "$RUNS_DIR"/longJobs_*.txt
		
		# Check the number of depth layers file exist
		if [[ ! -s "$INPUTDATA_DIR/$FILENAME_NUM_DEPTHLAYERS" ]]; then
    		die "$FILENAME_NUM_DEPTHLAYERS is missing or empty"
		fi

		# Read depth layers and categorise jobs
		jobId=1
		while IFS= read -r depthLayers; do
			if [ "$depthLayers" -le 200 ]; then
				shortJobsCount=$((shortJobsCount + 1))
				printf "%s\n" "$jobId" >> "$RUNS_DIR/shortJobs_${shortBatch}.txt"
				if [ "$shortJobsCount" -ge "$batchSizeSlurm" ]; then
					shortJobsCount=0
					shortBatch=$((shortBatch + 1))
				fi
			else
				longJobsCount=$((longJobsCount + 1))
				printf "%s\n" "$jobId" >> "$RUNS_DIR/longJobs_${longBatch}.txt"
				if [ "$longJobsCount" -ge "$batchSizeSlurm" ]; then
					longJobsCount=0
					longBatch=$((longBatch + 1))
				fi
			fi
			((jobId++))
		done < "$INPUTDATA_DIR/$FILENAME_NUM_DEPTHLAYERS"

		log "Short jobs split into $shortBatch batches."
		log "Long jobs split into $longBatch batches."
	fi
}

# --------------------------------------------------
# Submit job arrays
# --------------------------------------------------

submit_jobs_to_slurm_efficiently () 
{
    local job_type=$1   # short or long
    local mem=$2        # memory allocation for SLURM
    local time=$3       # time allocation for SLURM
    local partition=$4  # SLURM partition (short, long)
    local batch_file=$5 # path to the file containing job IDs to submit
    local batch_id=$6

	# Read job IDs into an array
    mapfile -t jobArray < "$batch_file"

    # Split job IDs into chunks of 1000
    runIdOffset=0
    batchSize=1000
    
    while [ "$runIdOffset" -lt "${#jobArray[@]}" ]; do

        # Get current batch of job IDs
        jobBatch=("${jobArray[@]:runIdOffset:batchSize}")
        lastArray="${#jobBatch[@]}"

        # Save job mapping to a file (for lookup inside jobs)
        mappingFile="job_mapping_${job_type}_${batch_id}_${runIdOffset}.txt"
		mappingPath="$RUNS_DIR/$mappingFile"
		rm -f "$mappingPath"
		: > "$mappingPath"
		for i in "${!jobBatch[@]}"; do
    		printf "%d %s\n" "$((i+1))" "${jobBatch[$i]}" >> "$mappingPath"
		done

        log "Submitting a $job_type job array of size $lastArray (batch offset: $runIdOffset)"
		unset SLURM_MEM_PER_NODE # ensure that --mem-per-cpu doesn't conflict with inherited memory settings
        sbatch --export=ALL,CONFIG_NAME="$CONFIG_NAME",OFFSET="$runIdOffset",MAPPING_FILE="$mappingPath",RUNS_DIR="$RUNS_DIR" \
               --time="$time" \
               --array=1-"$lastArray" \
               --mem-per-cpu="$mem" \
               --partition="$partition" \
               --qos=earth \
               --job-name=runscript \
               --output="$LOGS_DIR/log_runscript_%A_%a.log" \
               "$PIPELINE_DIR/runscript_balanced.sh"

        # Move to the next batch
        runIdOffset=$((runIdOffset + batchSize))

    done
}

# --------------------------------------------------
# Execution manager
# --------------------------------------------------

manage_runs ()
{
	# SLURM
	if [[ "$CHOICE_RUN_IN_SLURM" == "true" ]]; then

		# Short Jobs (≤5h)
		batch_counter=0
        for file in "$RUNS_DIR"/shortJobs_*.txt; do
            [[ -s "$file" ]] || continue  # skip if no files match
            submit_jobs_to_slurm_efficiently "short" "500MB" "03:00:00" "short" "$file" "$batch_counter"
            batch_counter=$((batch_counter+1))
        done

        # Long Jobs (>5h)
        batch_counter=0
        for file in "$RUNS_DIR"/longJobs_*.txt; do
            [[ -s "$file" ]] || continue  # skip if no files match
            submit_jobs_to_slurm_efficiently "long" "1000MB" "06:00:00" "long" "$file" "$batch_counter"
            batch_counter=$((batch_counter+1))
        done
        
		log "All SLURM job arrays submitted"
	
	# Local machine - a few locations
	elif [ "$CHOICE_RUN_IN_SLURM" == "false" ] && [ "$CHOICE_GRID_DOMAIN" == "2" ]; then
		
		for (( r=1; r<=nLocDirectories; r++ )); do		
			RUN_ID="$RUNS_DIR/run_$r"
			thisRunDir="$RUNS_DIR/run_$RUN_ID"
			cd "$thisRunDir"; pwd
			./SLAMSexecutable > log_"$RUN_ID" &
		done
	
	else
    	log "Invalid combination of CHOICE_RUN_IN_SLURM and CHOICE_GRID_DOMAIN"
	fi
}

# --------------------------------------------------
# Execution pipeline
# --------------------------------------------------
# Defines ordered workflow steps.

main () 
{
	log "Starting run_model_balanced pipeline"
	validate_config_name
	slams_paths_config
    slams_workflow_config
    slams_namelist_to_shell "$CONFIG_NAME" "WorkflowOptions" "WorkflowFilenames"
	count_number_location_directories
	determine_timing_of_jobs_in_slurm
	manage_runs
	log "Finished run_model_balanced pipeline"
}

# Execute this script
main