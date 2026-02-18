#!/bin/bash

# ======================================================================================
# Author: A. Rufas | 14 Feb 2026
#
# This script download runsoutput.mat files afetr sensitivity tests. It is to be run locally
#
# USAGE
# -----
#   ./run_param_sensitivity_test.sh
# ======================================================================================

# Define remote user and base path
USER="wolf4894"
SERVER="htc-login.arc.ox.ac.uk"
REMOTE_BASE="/data/eart-slam-dunk/wolf4894/Projects/SLAMS2.0/tests"
LOCAL_BASE="/Users/Anna/LocalDocuments/Academic/Projects/SLAMS2.0/tests"

# Get list of directories starting with LOCALTS6_ on the remote server
DIRS=$(ssh ${USER}@${SERVER} "ls -d ${REMOTE_BASE}/LOCALTS6_*/modelruns 2>/dev/null")

# Loop over each directory and copy the file
for REMOTE_MODEL_RUN in $DIRS; do
    # Extract the LOCALTS6_* directory name from the path
    DIR_NAME=$(basename $(dirname $REMOTE_MODEL_RUN))
    
    # Build local destination directory
    LOCAL_DEST="${LOCAL_BASE}/${DIR_NAME}/modelruns"
    
    # Create the local directory if it doesn't exist
    mkdir -p "$LOCAL_DEST"
    
    # Use scp to copy runsoutput.mat
    scp "${USER}@${SERVER}:${REMOTE_MODEL_RUN}/runsoutput.mat" "$LOCAL_DEST/"
done