#!/bin/bash

#SBATCH --job-name=cleanruns 
#SBATCH --time=10:00:00             
#SBATCH --partition=long			   
#SBATCH --mem-per-cpu=10GB   
#SBATCH --output=../../../logs/log_cleanruns_%A.log

set -euo pipefail # abort the script if errors with commands, variables and pipelines occur
IFS=$'\n\t'       # controls bash word splitting

# ====================================================================================== 
# SLAMS 2.0 – Input Preparation and Build Pipeline
#
# Author: A. Rufas | 16 Feb 2026    
#                                                                                                                                                                           
# PURPOSE
# -------
# Script to clean model run folders after detecting errors in model configuration files/code
# files and re-run if necessary.
#
# USAGE
# -----
#   Local:
#       ./clean_run_directories.sh <config_name> 
#
#   SLURM:
#       sbatch --export=ALL,CONFIG_NAME=<config_name> clean_run_directories.sh
#                                                                           
# ====================================================================================== 

# --------------------------------------------------
# Load core workflow modules
# --------------------------------------------------

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

# --------------------------------------------------
# Select fucntion(s) that need to be run
# --------------------------------------------------

main () 
{
	#slams_paths_config
    #slams_workflow_config
	remove_single_directory
	#remove_test_directories
	#remove_multiple_sequential_run_directories
	#remove_selected_run_directories
	#count_number_location_directories
	#remove_run_files
}

# --------------------------------------------------
# Functions used in this script
# --------------------------------------------------

# Remove a specific test directory and its subfolders
remove_single_directory ()
{
	log "Start cleaning"
	rm -rf /data/eart-slam-dunk/wolf4894/Projects/SLAMS2.0/tests/LOCALTS6_*
	log "Finished"
}

# Remove a sequence test directories and their subfolders
remove_test_directories () 
{
    for i in {100..161}; do
        rm -rf "/data/eart-slam-dunk/wolf4894/Projects/SLAMS2.0/tests/test$i"
    done
}

# Remove a range of run directories by iteration
remove_multiple_sequential_run_directories ()
{
	for (( i=4500; i<=5090; i++ )); do
		rundir="$RUNS_DIR/run_$i"
		echo "$rundir"
		rm -rf "$rundir" \;
	done
}

# Remove selected run directories
remove_selected_run_directories () 
{
    for i in 1 5; do
        rundir="$RUNS_DIR/run_$i"
        echo "Removing $rundir"
        rm -rf "$rundir"
    done
}

# Determine the number of water-column instances (or ocean locations) 
count_number_location_directories ()
{
	nLocDirectories=$(find "$INPUTDATA_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
	# Check if we have directories to process
    if [ "$nLocDirectories" -eq 0 ]; then
        echo "Error: No location-specific directories found."
        exit 1
    else
		echo "We have $nLocDirectories location-specific directories"
	fi
}

# Remove files that are used to run SLAMS in a run directory
remove_run_files ()
{
	echo "Removing run files..."	
	for (( i=1; i<="$nLocDirectories"; i++ )); do 
	#for i in 1 5; do
		rundir="$RUNS_DIR/run_$i"
		cd "$rundir"; pwd
		find . -type f -name "out_*" -exec rm -f {} \;
		find . -type f -name "log_*" -exec rm -f {} \;
		#find . -type f -name "namelist*" -exec rm -f {} \;
		#find . -type f -name "*.F90" -exec rm -f {} \;
		#find . -type f -name "*.mod" -exec rm -f {} \;
		#find . -type f -name "*.o" -exec rm -f {} \;
		#find . -type f -name "*.h" -exec rm -f {} \;
		#find . -type f -name "*.F" -exec rm -f {} \;
		#find . -type f -name "Makefile" -exec rm -f {} \;
		#find . -type f -name "SLAMS*" -exec rm -f {} \;
	done; echo "...all run directories cleaned."	
}

# Execute this script
main