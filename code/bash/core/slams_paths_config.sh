#!/usr/bin/env bash

# ====================================================================================== 
# SLAMS 2.0 – Define Project Directory Structure
#
# Author: A. Rufas | 16 Feb 2026    
#                                                                                                                                                                           
# PURPOSE
#   Define directory paths derived from ROOT_DIR and CONFIG_NAME.
#
# NOTE
#   This file defines reusable functions for SLAMS workflows.
#   It is not intended to be executed directly.
#   It is sourced by higher-level workflow scripts.
#   
# ======================================================================================                

slams_paths_config ()
{
	if [[ -z "${CONFIG_NAME:-}" ]]; then
        echo "CONFIG_NAME must be defined before calling define_project_paths" >&2
        return 1
    fi
    
	ROOT_DIR=/data/eart-slam-dunk/wolf4894/Projects/SLAMS2.0
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