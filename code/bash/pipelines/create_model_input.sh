#!/usr/bin/env bash
#SBATCH --job-name=input 
#SBATCH --time=06:00:00                   
#SBATCH --partition=short
#SBATCH --qos=standard
#SBATCH --mem=100GB
#SBATCH --output=../../../logs/log_input_%A.log

set -euo pipefail # abort the script if errors with commands, variables and pipelines occur
IFS=$'\n\t'       # controls bash word splitting

# ====================================================================================== 
# SLAMS – Input Preparation and Build Pipeline
#
# Author: A. Rufas
# Created: 16 Feb 2026    
# Last updated: 14 Sep 2026   
#                                                                                                                                                                           
# PURPOSE
# -------
#   Prepare a complete SLAMS experiment directory:
#     - Set up directory structure
#     - Compile model executable
#     - Generate or import forcing data
#     - Create per-run namelists
#
# USAGE
# -----
#   Local:
#       ./create_model_input.sh <config_name>
#
#   SLURM:
#       sbatch --export=ALL,CONFIG_NAME=<config_name> create_model_input.sh
#
# EXPECTED RUNTIME (approximate)
# ------------------------------
#   Global configuration                                : ~6 hours
#   Local configuration generating forcing data de novo : ~40 minutes
#   Local configuration reusing forcing data            : ~10 minutes
#
# NOTES
# -----
#   - Forcing data generation is the most computationally expensive step.
# 	- To avoid regeneration, set FORCING_MODE to import or reuse in slams_workflow_config.sh.
#
# ARCHITECTURE PRINCIPLE
# ----------------------
#   - This script defines execution workflow only.
#   - Scientific configuration parameters belong in the namelist, not here.
#   - The Makefile defines compilation behaviour.
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

source "$CORE_DIR/slams_logging.sh"
source "$CORE_DIR/slams_paths_config.sh"
source "$CORE_DIR/slams_workflow_config.sh"
source "$CORE_DIR/slams_environment_init.sh"
source "$CORE_DIR/slams_namelist_to_shell.sh"

# --------------------------------------------------
# Validate configuration terms
# --------------------------------------------------

validate_configuration ()
{
	log "About to validate configuration"
	
    [[ -n "${FORCING_MODE:-}" ]] \
        || die "FORCING_MODE is not set"

    case "$FORCING_MODE" in
        import)
            [[ -n "${IMPORT_FROM_TESTDIR:-}" ]] \
                || die "IMPORT_FROM_TESTDIR must be set when FORCING_MODE=import"

            [[ "$IMPORT_FROM_TESTDIR" == "$CONFIG_NAME" ]] \
                && die "Cannot import input data from the same test directory"
            ;;
        regenerate|reuse)
            ;;
        *)
            die "Invalid FORCING_MODE: $FORCING_MODE"
            ;;
    esac
    
    log "Validation completed"
}

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
# Prepare the directory structure for model runs
# --------------------------------------------------
# Ensures that:
#   - Required directory structure exists
#   - Old run artifacts are removed (logs and run directories)

setup_workspace () 
{
	mkdir -p "$TEST_DIR"
    mkdir -p "$TEST_DIR/logs"
    mkdir -p "$TEST_DIR/modelruns"
    mkdir -p "$MODEL_DIR"

    rm -rf "$TEST_DIR/logs/"*
    rm -rf "$TEST_DIR/modelruns/"*

    mkdir -p "$INPUTDATA_DIR"
}

# --------------------------------------------------
# Compile SLAMS executable
# --------------------------------------------------
# Optionally re-syncs source tree into build directory and compiles using Makefile.
# Compiler selection depends on CHOICE_RUN_IN_SLURM.

compile_code ()
{
	if [ "$CHOICE_RECOMPILE" == "true" ]; then
		log "Source code has been modified --> syncing source tree"
    	rm -rf "$MODEL_DIR"/*
    	cp -r "$SOURCE_DIR"/* "$MODEL_DIR"

		# Go to model dir
		pushd "$MODEL_DIR" >/dev/null || die "Cannot enter $MODEL_DIR"
		log "Compiling model code in $(pwd)"
	
		# Ensure a clean compile
    	make clean || true
    	make cleanoutput || true
    
		# Compile
    	if [[ "$CHOICE_RUN_IN_SLURM" == "true" ]]; then
    		# Check the root directory of the MKL installation managed by EasyBuild
        	[[ -n "${EBROOTIMKL:-}" ]] || {
            	die "EBROOTIMKL not set (imkl module not loaded)"
        	}
        	# Check the directory that contains MKL .so files exists (MKL_LIBDIR is initialised in slams_preamble)
        	[[ -n "${MKL_LIBDIR:-}" ]] || {
    			die "MKL_LIBDIR not set (MKL not correctly initialised)"
			}
        	MKL_RPATH="-Wl,-rpath,${MKL_LIBDIR}"
        	make FC=ifort DEBUG="$DEBUG" MKL_RPATH="$MKL_RPATH" -j1
    	else
        	if command -v xcrun >/dev/null 2>&1; then
            	export SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"
        	fi
        	make FC=gfortran DEBUG="$DEBUG" -j1
    	fi

		# Return to previous directory
    	log "Compilation done"
    	popd >/dev/null    
    fi
}

# --------------------------------------------------
# Prepare forcing input directory
# --------------------------------------------------
# Cleans or preserves INPUTDATA_DIR according to FORCING_MODE.
# Separation from setup_workspace() avoids coupling forcing logic with infrastructure setup.

prepare_inputdata_directory () 
{
    if [[ ! -d "$INPUTDATA_DIR" ]]; then
        mkdir -p "$INPUTDATA_DIR"
        return
    fi

    log "Preparing input data directory"

    case "$FORCING_MODE" in
        regenerate)
            find "$INPUTDATA_DIR" -type f -name "$FILENAME_NUM_DEPTHLAYERS" -delete
            find "$INPUTDATA_DIR" -type f -name "$FILENAME_RUN_GRID" -delete
            find "$INPUTDATA_DIR" -type f -name "lonlatRef.txt" -delete
            find "$INPUTDATA_DIR" -type f -name "namelist*" -delete
            find "$INPUTDATA_DIR" -depth -type d -name "run*" -exec rm -rf {} +
            ;;
        import)
            # Only remove transient files; do not remove run_* yet
            find "$INPUTDATA_DIR" -type f -name "$FILENAME_RUN_GRID" -delete
            find "$INPUTDATA_DIR" -type f -name "lonlatRef.txt" -delete
            find "$INPUTDATA_DIR" -type f -name "namelist*" -delete
            ;;
        reuse)
            log "Keeping existing input data directory unchanged"
            ;;
    esac
}

# --------------------------------------------------
# Grid selection (MATLAB)
# --------------------------------------------------
# Generates run grid from template according to CHOICE_GRID_DOMAIN.

select_grid ()
{
	log "Creating grid using the template $FILENAME_TEMPLATE_GRID"
	matlab -nodesktop -nosplash -r \
		"addpath(genpath('$MATLABCODE_DIR/')); \
		selectGridDomain('$RAW_DATA_DIR/','$INPUTDATA_DIR/','$CONFIG_DIR/','$FILENAME_TEMPLATE_GRID',\
						 '$FILENAME_RUN_GRID','$FILENAME_LOCAL_LATITUDES','$FILENAME_LOCAL_LONGITUDES',\
						 $CHOICE_GRID_DOMAIN,[]); \
		clc; clear all; exit(0);" || die "MATLAB execution failed"
	
	# Now verify expected output exists
	if [[ ! -f "$INPUTDATA_DIR/$FILENAME_RUN_GRID" ]]; then
    	die "Grid file not generated: $INPUTDATA_DIR/$FILENAME_RUN_GRID"
	fi
}    

# --------------------------------------------------
# Forcing lifecycle dispatcher
# --------------------------------------------------
# Executes forcing workflow according to FORCING_MODE.
# Ensures that only one forcing path is taken.

handle_forcing_strategy () 
{
    case "${FORCING_MODE}" in
        import)
        	prepare_inputdata_directory
            log "Importing forcing data from ${IMPORT_FROM_TESTDIR}"
            import_model_forcing_files
            ;;
        regenerate)
            log "Preparing input directory for regeneration"
            prepare_inputdata_directory
            log "Creating grid"
            select_grid
            log "Generating new forcing input data"
            create_model_forcing_files
            ;;
        reuse)
            log "Reusing existing local forcing input data"
            ;;
        *)
            die "Unknown FORCING_MODE: ${FORCING_MODE}"
            ;;
    esac
}

# --------------------------------------------------
# Import pre-generated forcing data (optional)
# --------------------------------------------------
# Copies previously generated forcing files into the current experiment directory to 
# avoid rerunning the expensive MATLAB preprocessing step.

import_model_forcing_files ()
{
    SRC_INPUTDATA="$ROOT_DIR/tests/${IMPORT_FROM_TESTDIR}/modelinputdata"

    if [[ ! -d "$SRC_INPUTDATA" ]]; then
        die "Source input data directory not found: $SRC_INPUTDATA"
    fi

    # Create target input data directory if needed
    mkdir -p "$INPUTDATA_DIR"

    # Copy only what matters
    rsync -av \
        --include="run_*/" \
        --include="run_*/*" \
        --include="$FILENAME_NUM_DEPTHLAYERS" \
        --include="$FILENAME_RUN_GRID" \
    	--include="slamsForcing.mat" \
    	--include="lonlatRef.txt" \
        --exclude="*" \
        "$SRC_INPUTDATA/" "$INPUTDATA_DIR/"

    log "Forcing data imported successfully"
}

# --------------------------------------------------
# Generate forcing data (MATLAB)
# --------------------------------------------------
# Runs createSlamsInputData.m to produce binary input files required by SLAMS executable.

create_model_forcing_files ()
{
	log "Creating model input forcing files and sending them to their corresponding directories"
	matlab -nodesktop -nosplash -r \
		"addpath(genpath('$MATLABCODE_DIR/')); \
		addpath('$RESOURCES_INTERNAL_DIR/'); \
		addpath(genpath('$RESOURCES_EXTERNAL_DIR/')); \
		addpath('$RAW_DATA_DIR/'); \
		addpath('$INPUTDATA_DIR/'); \
		addpath('$INTERIM_DATA_DIR/'); \
		createSlamsInputData('$RAW_DATA_DIR/','$INPUTDATA_DIR/','$INTERIM_DATA_DIR/','$FILENAME_DATA_PAR0',\
							 '$FILENAME_DATA_NPP','$FILENAME_DATA_CHLA','$FILENAME_DATA_RHO','$FILENAME_DATA_OMEGACALC',\
							 '$FILENAME_DATA_NIT','$FILENAME_DATA_SIL','$FILENAME_DATA_PHOS','$FILENAME_DATA_OXY',\
							 '$FILENAME_DATA_TEMP','$FILENAME_DATA_DUST','$FILENAME_DATA_MLD','$FILENAME_DATA_MESOZOO',\
							 '$FILENAME_DATA_DYNVISCO','$FILENAME_DATA_NUM_DAYLIGHTHOURS','$FILENAME_DATA_MASK',\
							 '$FILENAME_RUN_GRID','$FILENAME_NUM_DEPTHLAYERS',$CHOICE_INPUT_DATA_TO_FORCE_MODEL,\
							 $CHOICE_LIGHT_DISTRIBUTION,$CHOICE_GRID_DOMAIN); \
		clc; clear all; exit(0);" || die "MATLAB execution failed"
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
# Transfer code and input files
# --------------------------------------------------
# For each run_i directory:
#   - Create run directory
#   - Link input forcing data
#   - Copy executable and namelist

transfer_code_and_input_files ()
{	
	log "Transferring code files and input data files"

	# Iterate over location-specific directories
	for (( i=1; i<=nLocDirectories; i++ )); do
		currdir="$RUNS_DIR/run_$i"
		mkdir "$currdir"

		# Soft link to input data files
		ln -s "$INPUTDATA_DIR/run_$i"/* "$currdir" 
		
		# Hard link to model files (executable)
		cp "$MODEL_DIR/SLAMSexecutable" "$currdir/"
		
		# Write processed namelist directly
  		printf "%s\n" "$NAMELIST_TEMP" > "$currdir/namelist.input"
	done
	
	log "Done"
}

# --------------------------------------------------
# Adjust per-run namelist parameters
# --------------------------------------------------
# Updates nDepthLayers in each run directory according to forcing-specific vertical resolution.

create_input_namelist ()
{
    # Ensure the input file exists
    if [ ! -f "$INPUTDATA_DIR/$FILENAME_NUM_DEPTHLAYERS" ]; then
        die "File '$INPUTDATA_DIR/$FILENAME_NUM_DEPTHLAYERS' not found"
    fi
    
	# Read in the entries in FILENAME_NUM_DEPTHLAYERS (created by createSlamsInputData.m)
	dataNumDepthLayers=() # create an empty array
	while IFS= read -r line; do # reading file in row mode, insert each line into array
		dataNumDepthLayers+=("$line") # append
	done < "$INPUTDATA_DIR/$FILENAME_NUM_DEPTHLAYERS"

    # Check if the number of elements in dataNumDepthLayers matches the number of directories
    if [ "${#dataNumDepthLayers[@]}" -ne "$nLocDirectories" ]; then
        die "The number of depth layer entries does not match the number of directories"
    fi
    
	# Loop through the run directories and replace the default nDepthLayers with the 
	# actual value for that run instance (now, in dataNumDepthLayers)
	log "Writing the local number of depth layers into the namelist"
	i=0 # 'i' is for going through dataNumDepthLayers array
	for (( j=1; j<=nLocDirectories; j++ )); do
		currdir="$RUNS_DIR/run_$j"

		# Prepare the replacement string
		replacementString="nDepthLayers=${dataNumDepthLayers[$i]}"
		log "The replacement string is $replacementString"
		
		# Replace the string in the namelist input
		sed -i "s/nDepthLayers=.*/$replacementString/" "$currdir/namelist.input"
		((i=i+1))
	done
	log "Done"
}

# --------------------------------------------------
# Smoke test
# --------------------------------------------------
# Executes a short non-blocking run to verify that the executable launches and enters main().

smoke_test ()
{
	set +e   # <-- CRITICAL: allow non-zero exits in smoke test
	
	local SAMPLE_RUN_DIR="${RUNS_DIR}/run_1"	
	log "Running smoke test in ${SAMPLE_RUN_DIR} (this should start the program and reveal runtime errors)"
	
	if [[ ! -d "${SAMPLE_RUN_DIR}" ]]; then
        die "Sample run dir ${SAMPLE_RUN_DIR} does not exist"
    fi
    
    pushd "${SAMPLE_RUN_DIR}" >/dev/null

    if [[ ! -x "./SLAMSexecutable" ]]; then
    	popd >/dev/null
        die "SLAMSexecutable missing or not executable in ${SAMPLE_RUN_DIR}"
    fi

 	log "---- smoke test start ----"
 	
 	./SLAMSexecutable </dev/null > smoke.out 2>&1 &
	SMOKE_PID=$!
	
	# wait up to 10 seconds
	for i in {1..10}; do
		if ! kill -0 $SMOKE_PID 2>/dev/null; then
			break
		fi
		sleep 1
	done
	# If still running after timeout → kill
    if kill -0 "$SMOKE_PID" 2>/dev/null; then
        kill -9 "$SMOKE_PID"
        wait "$SMOKE_PID" 2>/dev/null
        popd >/dev/null
        set -e
        die "Smoke test hung — process killed"
    fi

    wait "$SMOKE_PID"
    SMOKE_EXIT=$?

    SMOKE_OUT=$(cat smoke.out)

    log "---- smoke test end (exit=${SMOKE_EXIT}) ----"
    log "${SMOKE_OUT}"

    # Check that program reached main()
    if echo "$SMOKE_OUT" | grep -qi "The SLAMS program has been entered successfully"; then
        log "Smoke test PASSED: executable entered main()"
    else
        popd >/dev/null
        set -e
        die "Smoke test FAILED — executable did not reach main()"
    fi

    popd >/dev/null
    set -e
}

# --------------------------------------------------
# Execution pipeline
# --------------------------------------------------
# Defines ordered workflow steps.
# Modify PIPELINE array to enable/disable stages.

main () 
{
	log "Starting create_model_input pipeline"
	
	validate_config_name
	
	slams_paths_config
    slams_workflow_config
    
	slams_environment_init
	slams_namelist_to_shell "$CONFIG_NAME" "WorkflowOptions" "ForcingDataFiles" "WorkflowFilenames"

	validate_configuration
	
    PIPELINE=(
        setup_workspace
        compile_code
        handle_forcing_strategy
        count_number_location_directories
        transfer_code_and_input_files
        create_input_namelist
        # smoke_test   # enable if needed
    )
    
    for step in "${PIPELINE[@]}"; do
        log "Executing: $step"
        "$step"
    done
    
 	# Minimal debug summary 
 	log "CHOICE_GRID_DOMAIN: $CHOICE_GRID_DOMAIN"
    log "TEST_DIR: $TEST_DIR"
    log "MKL_LIBDIR=${MKL_LIBDIR:-<unset>}"
    log "INTEL_OMP_LIB=${INTEL_OMP_LIB:-<unset>}"
    log "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-<unset>}"

 	log "Workspace ready"
}

# Execute this script
main 	