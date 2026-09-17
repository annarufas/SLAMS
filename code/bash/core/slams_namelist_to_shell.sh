#!/usr/bin/env bash

# ====================================================================================== 
# SLAMS – Namelist Section Parser
#
# Author: A. Rufas 
# Created: 16 Feb 2026    
# Last updated: 14 Sep 2026
#                                                                                                                                                                           
# PURPOSE
# -------
#   Parse selected sections of a SLAMS namelist file and export the requested variables 
#   into the current shell environment.
#
# NOTES
# -----
# 	- This file contains internal functions used by SLAMS workflows.
#	- It is sourced by higher-level pipeline scripts and is not intended to be executed or
#   modified directly by users.
# 
# ======================================================================================                

slams_namelist_to_shell () 
{
    local CONFIG_NAME=$1
    shift
    local NAMELIST_SECTIONS_TO_READ=("$@")  # array of namelist section names to read

    # ------- Locate namelist -------
    local NAMELIST_FILE="${ROOT_DIR}/config/namelist_${CONFIG_NAME}.txt"
    if [[ ! -f "$NAMELIST_FILE" ]]; then
        die "Namelist file not found: $NAMELIST_FILE"
    fi
    log "Parsing namelist: $NAMELIST_FILE"
    
    # ---- Process namelist in memory ----
	NAMELIST_TEMP=$(envsubst < "$NAMELIST_FILE" \
        | sed 's/^[ \t]*//;s/[ \t]*$//') \
        || die "Failed to process namelist"

	# Guard against unresolved variables
    if grep -q '\$[A-Za-z_][A-Za-z0-9_]*' <<< "$NAMELIST_TEMP"; then
        die "Unresolved variables remain in processed namelist"
    fi
    
	# ------- Parse requested sections -------
    local current_section=""
    local in_target_section=false
    
    while IFS= read -r line || [[ -n "$line" ]]; do # read the namelist line by line
    
        # Check for the start of a section (e.g., &General)
        if [[ "$line" =~ ^\&([A-Za-z0-9_]+) ]]; then
            current_section="${BASH_REMATCH[1]}"
            in_target_section=false

            for section in "${NAMELIST_SECTIONS_TO_READ[@]}"; do
                if [[ "$current_section" == "$section" ]]; then
                    in_target_section=true
                    break
                fi
            done
        fi
    
        # Check for the end of a section (indicated by a '/')
        if [[ "$line" =~ ^/ ]]; then
            in_target_section=false
        fi
    
        # Process lines within the target sections
        if $in_target_section; then
            # Skip comments and empty lines
            [[ "$line" =~ ^[[:space:]]*! ]] && continue
            [[ -z "$line" ]] && continue
    
            # Match variable definitions and export them
            if [[ "$line" =~ ^([A-Za-z0-9_]+)[[:space:]]*=[[:space:]]*(.+)[[:space:]]*!?.*$ ]]; then
                local key="${BASH_REMATCH[1]}"
                local value="${BASH_REMATCH[2]}"
    
                # Remove leading and trailing spaces
                value=$(echo "$value" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    
                # Handle special cases for `.true.` and `.false.`
                if [[ "$value" == ".true." ]]; then
                    value=true
                elif [[ "$value" == ".false." ]]; then
                    value=false
                fi
    
                # Export the variable
                export "$key=$value"
            fi
        fi
    done <<< "$NAMELIST_TEMP"
    log "Namelist sections parsed successfully"
}
