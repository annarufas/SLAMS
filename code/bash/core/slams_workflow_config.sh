#!/usr/bin/env bash

# ====================================================================================== 
# SLAMS – Define Workflow Options
#
# Author: A. Rufas
# Created: 16 Feb 2026    
# Last updated: 14 Sep 2026  
#                                                                                                                                                                           
# PURPOSE
# -------
#   Defines workflow options that control how the SLAMS pipeline operates.
#
# NOTES
# -----
#   - This file defines the workflow configuration function used by SLAMS pipelines.
#   - It is sourced by pipeline scripts, not executed directly. 
#
# USER CONFIGURATION
# ------------------
# 	Edit the default values below to configure a particular SLAMS experiment.
#   
# ======================================================================================                

#   CHOICE_RUN_IN_SLURM : execution mode
#							true  = run under SLURM
#							false = run on a local machine
#   CHOICE_RECOMPILE    : force recompilation 
#							true  = force recompilation
#							false = reuse the existing executable when available
#   FORCING_MODE        : forcing-data strategy
#							regenerate = rebuild forcing data from raw data
#							import     = copy forcing from another experiment
#							reuse      = use existing local forcing data
#   IMPORT_FROM_TESTDIR : source experiment directory when FORCING_MODE=import
#   DEBUG               : Makefile build mode
#							1 = debug mode active
#							0 = debug mode inactive

slams_workflow_config ()
{
	CHOICE_RUN_IN_SLURM=true
	CHOICE_RECOMPILE=true
	FORCING_MODE=regenerate
	IMPORT_FROM_TESTDIR=test167
	DEBUG="${DEBUG:-0}"
}