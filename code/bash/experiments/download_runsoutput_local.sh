#!/usr/bin/env bash

# ======================================================================================
# SLAMS – Download Sensitivity-Test Outputs
#
# Author: A. Rufas
# Created: 16 Feb 2026
# Last updated: 14 Sep 2026
#
# PURPOSE
# -------
# 	Downloads runsoutput.mat files from sensitivity-test runs on the ARC HPC system to the
# 	corresponding local SLAMS test directories.
#
# USAGE
# -----
# 	./download_runsoutput_local.sh
#
# REQUIREMENTS
# ------------
# 	- SSH access to the ARC HPC system
# 	- scp
#
# USER CONFIGURATION
# ------------------
# 	Edit the remote and local paths below if your directory structure differs.
#
# ======================================================================================

# Define remote user and base path
USER="wolf4894"
SERVER="htc-login.arc.ox.ac.uk"
REMOTE_BASE="/data/eart-slam-dunk/wolf4894/Projects/SLAMS/tests"
LOCAL_BASE="/Users/Anna/LocalDocuments/Academic/Projects/SLAMS/tests"

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