# ........................................................................................
#
# Makefile for compiling and linking SLAMS-2.0 (Fortran90).
#
# NOTE:
# This Makefile assumes that the required compiler toolchain and external libraries
# (e.g., Intel compiler, MKL) are already available in the environment.
#
# On HPC systems, this is typically handled via module load commands in a SLURM
# submission script. In this project, environment setup is performed externally
# (e.g. in create_model_input.sh) before invoking make.
#
# The Makefile itself does not load modules or configure library paths.
#
# The build system relies on a manually maintained list of source files. Failure to update 
# this list when modifying the source tree may result in compilation errors or inconsistent 
# builds.
#
# ........................................................................................

# ======================================================
# Compiler selection
# ======================================================

# Fortran compiler (can be overridden at command line)
# Example:
#   make FC=gfortran ! GNU's Fortran90 compiler
#   make FC=ifort    ! Intel's Fortran90 compiler

FC ?= gfortran 

# ======================================================
# Linker configuration
# ======================================================

# Additional linker flags (e.g. library search paths)
# Example:
#   make LDFLAGS="-L/usr/local/lib"

LDFLAGS ?=

# External libraries to link against
# Example:
#   make LIBS="-lm -lblas -llapack"

LIBS ?=

# Runtime library search path (mainly used with Intel MKL)
# Example:
#   make MKL_RPATH="-Wl,-rpath,/opt/intel/mkl/lib/intel64"

MKL_RPATH ?=

# ======================================================
# Build configuration
# ======================================================

# DEBUG mode:
#   DEBUG=1  ! development mode (no optimisation + runtime checks)
#   DEBUG=0  ! production mode (optimised build)

DEBUG ?= 0

# Optimisation flags
ifeq ($(DEBUG),1)
	OPTIM_FLAGS = -g -O0
else
	OPTIM_FLAGS = -g -O2 # -O3 higher level of optimisation compared to -O2
endif

# ======================================================
# GNU Fortran (gfortran) configuration
# ======================================================

ifeq ($(FC),gfortran)
	
	# Compiler warning flags (always enabled)
	WARN_FLAGS = -Wall -Wextra -Wconversion

	# Additional runtime checking in DEBUG mode
  	ifeq ($(DEBUG),1)
    	CHECK_FLAGS = -fcheck=all -Wuninitialized -Wunused -Wsurprising -Waliasing \
                 	  -ffpe-trap=invalid,zero,overflow -fimplicit-none -fbacktrace
  	else
    	CHECK_FLAGS =
  	endif
	
	# Final compiler flags
  	FCFLAGS = $(OPTIM_FLAGS) $(WARN_FLAGS) $(CHECK_FLAGS) -fconvert=big-endian
 	MOD_DIR_FLAG = -J .
 	
endif

# ======================================================
# Intel Fortran (ifort) configuration
# ======================================================

ifeq ($(FC),ifort)
	
	# Enable broad compile-time diagnostics
  	WARN_FLAGS = -warn all
	
	# Additional runtime checking in DEBUG mode
  	ifeq ($(DEBUG),1)
    	CHECK_FLAGS = -check all -check bounds -check uninit -warn interfaces -fpe0 -ftrapuv \
                      -init=snan -gen-interfaces
  	else
    	CHECK_FLAGS =
  	endif

	# Final compiler flags
  	FCFLAGS = $(OPTIM_FLAGS) $(WARN_FLAGS) $(CHECK_FLAGS) -fpp -convert big_endian \
              -assume byterecl -traceback -qmkl
  	MOD_DIR_FLAG = -module .
  	
endif

# ======================================================
# Source file definitions
# ======================================================

# List of Fortran 90 source files.
# 
# IMPORTANT:
# If any .F90 is renamed, removed, or newly added, the corresponding entry must be updated 
# in the SRCF90 variables below.
#
# The order of files matters: modules must appear before any files that use them.

SRCF90 = particlestructure.F90 \
		 safemath.F90 \
		 modelcounters.F90 \
		 modelconstants.F90 \
		 modelparameters.F90 \
		 modelgrid.F90 \
		 modelforcingdata.F90 \
		 modeleulerianvariables.F90 \
		 modelparticlecollection.F90 \
		 sanitychecks.F90 \
		 initialisation.F90 \
		 montecarlosampling.F90 \
		 findfunctions.F90 \
		 waterphysicsandlight.F90 \
		 calcparticleattributes.F90 \
		 mineraldissolution.F90 \
		 primaryproduction.F90 \
		 injectsurfaceparticles.F90 \
		 particledynamics.F90 \
		 heterotrophicmetabolism.F90 \
		 sink.F90 \
		 modeloutput.F90 \
		 wrappers.F90 \
		 timer.F90 \
		 theseed.F90 \
		 main.F90
SRCF = io.F

# ======================================================
# Derived build objects
# ======================================================
		
# Automatically generate object file names from source files
# (.F90 → .o, .F → .o)

OBJS = $(SRCF90:.F90=.o) $(SRCF:.F=.o)

# Name of final executable

EXEC = SLAMSexecutable

# Header file containing C-style preprocessor definitions (included via #include statements)

BDF = blockdefinitions.h

# ======================================================
# Build targets
# ======================================================

# Links all compiled object files into a single executable.

$(EXEC): $(OBJS) $(BDF)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(MKL_RPATH) -o $(EXEC) $(OBJS) $(LIBS)

# ======================================================
# Compilation rules
# ======================================================

%.o : %.F90 $(BDF)
	$(FC) $(FCFLAGS) $(MOD_DIR_FLAG) -c $< -o $@

%.o : %.F $(BDF)
	$(FC) $(FCFLAGS) $(MOD_DIR_FLAG) -c $< -o $@
	
# ======================================================
# Cleaning targets
# ======================================================

.PHONY: clean

# Clean all object files, executable, and module files

clean:
	rm -f $(OBJS) $(EXEC) *.mod

# Clean output files (logs and output data)

cleanoutput:
	rm -f *.out log* out*
	