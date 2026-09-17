#!/usr/bin/env bash
#SBATCH --job-name=checkout 
#SBATCH --time=00:10:00                   
#SBATCH --partition=devel
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../../../logs/log_checkout_%A.log

set -euo pipefail # abort the script if errors with commands, variables and pipelines occur
IFS=$'\n\t'       # controls bash word splitting

# ====================================================================================== 
# SLAMS – Check Run Status
#
# Author: A. Rufas
# Created: 16 Feb 2026
# Last updated: 16 Sep 2026  
#                                                                                                                                                                           
# PURPOSE
# -------
#   Inspects completed run directories and resubmits incomplete or failed runs when needed.
#
# USAGE
# -----
#   Local:
#       ./check_output_status.sh <config_name> 
#
#   SLURM:
#       sbatch --export=ALL,CONFIG_NAME=<config_name> check_output_status.sh
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
source "$CORE_DIR/slams_environment_init.sh"

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
# Inspect each run directory and decide the action
# --------------------------------------------------

check_status_run ()
{
	set +e  # disable automatic exit on error
	nFile=0

	for (( i=1; i<="$nLocDirectories"; i++ )); do 
		rundir="$RUNS_DIR/run_$i"
		controlfile="$rundir/out_control.bin"
    	statusfile="$rundir/program_status.tmp"
		
		if [[ ! -f "$controlfile" ]]; then # there's no control file
		
			if [[ -f "$statusfile" ]]; then
			
                status=$(cat "$statusfile")

                if [[ "$status" == "no clusters seeded" ]]; then
                    log "No clusters were seeded for run_$i. Skipping re-run."
                    continue
                elif [[ "$status" == "water too shallow" ]]; then
                    log "Water was too shallow for run_$i. Skipping re-run."
                    continue
                elif [[ "$status" == "in progress" ]]; then
                    log "Run_$i did not have time to finish. Resubmitting job with updated max. time..."
					submit_job "$i"
                    ((nFile++))
                    continue
                else
                    log "Unknown status in run_$i: $status. Investigating..."  
                    submit_job "$i"           
                    ((nFile++))
                    continue
                fi
                
            else
            
                log "No status nor control file exist in run_$i, indicating a problem with the job."
                submit_job "$i"
                ((nFile++))

            fi # checking if status file exists

		else # there's control file
		
			if [[ -f "$statusfile" ]]; then
                status=$(cat "$statusfile")
		
		        if [[ "$status" == "completed" ]]; then
                	log "Run_$i completed successfully."
                	continue
                else
                    log "Run_$i status is '$status'. Further investigation needed."
                    submit_job "$i"
                    ((nFile++))
                    continue
                fi
            fi
		
		fi
	done # loop over runs
	
	set -e  # Re-enable automatic exit on error
	log "Number of runs that require to be re-run: $nFile"
}

# --------------------------------------------------
# Submit SLURM job for specific run
# --------------------------------------------------

submit_job ()
{
    local run_id="$1"
    sbatch --export=ALL,CONFIG_NAME="$CONFIG_NAME",RUN_ID="$run_id",RUNS_DIR="$RUNS_DIR" \
           --time=10:00:00 \
           --array=0 \
           --mem-per-cpu=1000MB \
           --partition=short \
           --qos=standard \
           --job-name=runscript \
           --output="$LOGS_DIR/log_runscript_%A_%a.log" \
           "$PIPELINE_DIR/runscript_balanced.sh"
}

# --------------------------------------------------
# Execution pipeline
# --------------------------------------------------
# Defines ordered workflow steps.

main () 
{
	log "Starting run status check"
	slams_paths_config
	slams_environment_init
	count_number_location_directories
	check_status_run
 	log "Run status check complete"
}

# Execute this script
main