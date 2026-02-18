#!/usr/bin/env bash

# ====================================================================================== 
# SLAMS 2.0 – Define Workflow Options
#
# Author: A. Rufas | 16 Feb 2026    
#                                                                                                                                                                           
# PURPOSE
#   Controls how the pipeline behaves.
#
# NOTE
#   This file defines reusable functions for SLAMS workflows.
#   It is not intended to be executed directly.
#   It is sourced by higher-level workflow scripts.
#   
# ======================================================================================                

#   CHOICE_RUN_IN_SLURM : run mode
#							true = run under SLURM
#							false = run on local machine
#   CHOICE_RECOMPILE    : force recompilation 
#							true
#							false
#   FORCING_MODE        : forcing data strategy
#							regenerate = rebuild forcing data from raw
#							import = copy forcing from another experiment
#							reuse = use existing local forcing
#   IMPORT_FROM_TESTDIR : source experiment directory when FORCING_MODE=import
#   DEBUG               : Makefile build mode
#							1 = debug active
#							0 = debug inactive

slams_workflow_config ()
{
	CHOICE_RUN_IN_SLURM=true
	CHOICE_RECOMPILE=false
	FORCING_MODE=regenerate
	IMPORT_FROM_TESTDIR=test167
	DEBUG="${DEBUG:-0}"
}