#!/bin/bash

set -e  

module load intel/2022a 

# Define RUN_ID based on SLURM_ARRAY_TASK_ID and the mapping file
RUN_ID=$((OFFSET + SLURM_ARRAY_TASK_ID)) # $SLURM_ARRAY_TASK_ID is given by --array=1-$lastArray
#RUN_ID=$((SLURM_ARRAY_TASK_ID)) # $SLURM_ARRAY_TASK_ID is given by --array=1-$lastArray

# Define the directory for this specific run
thisRunDir="$RUNS_DIR"/run_"$RUN_ID"

# Ensure the directory exists and change to it
cd "$thisRunDir" || { echo "Failed to change directory to $thisRunDir"; exit 1; }
echo "Running job ID $RUN_ID in directory: $(pwd)"
make cleanoutput

# Run the executable
./SLAMSexecutable > log_"$RUN_ID" & # if you want background execution, keep the '&'
wait # wait for background jobs to complete before exiting
