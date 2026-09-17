#!/usr/bin/env bash

# ====================================================================================== 
# SLAMS – Define Project Directory Structure
#
# Author: A. Rufas
# Created: 16 Feb 2026    
# Last updated: 14 Sep 2026   
#                                                                                                                                                                           
# PURPOSE
# -------
#   Define directory paths derived from ROOT_DIR and CONFIG_NAME.
#
# NOTES
# -----
# 	- This file defines the directory-configuration function used by SLAMS pipelines.
#   - It is sourced by pipeline scripts, not executed directly. 
#
# USER CONFIGURATION
# ------------------
# 	- Set ROOT_DIR to the root directory of the SLAMS project.
# 	- The default value below is configured for my HPC environment.
#   
# ======================================================================================                

slams_paths_config ()
{
	if [[ -z "${CONFIG_NAME:-}" ]]; then
        echo "CONFIG_NAME must be defined before calling slams_paths_config" >&2
        return 1
    fi
    
    # User-configurable project root
	ROOT_DIR=/data/eart-slam-dunk/wolf4894/Projects/SLAMS
	
	TEST_DIR=$ROOT_DIR/tests/$CONFIG_NAME
	SOURCE_DIR=$ROOT_DIR/code/fortran/src
	MODEL_DIR=$ROOT_DIR/code/fortran/build
	RAW_DATA_DIR=$ROOT_DIR/data/raw
	PROCESSED_DATA_DIR=$ROOT_DIR/data/processed
	INTERIM_DATA_DIR=$ROOT_DIR/data/interim
	RESOURCES_INTERNAL_DIR=$ROOT_DIR/modelresources/internal
	RESOURCES_EXTERNAL_DIR=$ROOT_DIR/modelresources/external
	MATLABCODE_DIR=$ROOT_DIR/code/matlab
	CONFIG_DIR=$ROOT_DIR/config
	INPUTDATA_DIR=$ROOT_DIR/tests/$CONFIG_NAME/modelinputdata
	RUNS_DIR=$ROOT_DIR/tests/$CONFIG_NAME/modelruns
	LOGS_DIR=$ROOT_DIR/tests/$CONFIG_NAME/logs
	OUT_DIR=$ROOT_DIR/tests/$CONFIG_NAME/modeloutput
	OBS_DIR=$ROOT_DIR/tests/$CONFIG_NAME/observations
	PLOTS_DIR=$ROOT_DIR/tests/$CONFIG_NAME/results/$CONFIG_NAME
}