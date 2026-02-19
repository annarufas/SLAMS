# SLAMS-2.0

SLAMS-2.0 (**Stochastic, Lagrangian Aggregate Model of Sinking particles, v2.0**) is a fast, modular Fortran 90 particle-tracking model designed to simulate biogenic marine particle attributes and their dynamics in the ocean's biological carbon pump. 

## Key Features

- Stochastic Lagrangian particle-resolving framework
- Super-Droplet Method ([Shima et al., 2009](https://doi.org/10.1002/qj.441)) for efficient resolution of particle interaction dynamics
- Physically grounded particle evolution via fractal scaling and Reynolds-number-dependent settling
- Explicit resolution of aggregation, fragmentation, and zooplankton-mediated repackaging
- Particle fluxes and size spectra emerge diagnostically from explicit depth crossings
- One-dimensional water-column framework scalable to global grids via independent column execution
- Modular Fortran 90 core with Bash workflow control and MATLAB analysis tools

## Architecture

SLAMS-2.0 adopts a layered architecture separating scientific computation from workflow control:
1. Model layer (**Fortran 90**)
   - Core SLAMS-2.0 simulation engine
   - Reads binary forcing files
   - Scientific parameters configurable via `namelist`
2. Workflow orchestration layer (**Bash**)
   - Environment initialisation
   - Compilation
   - Forcing management
   - SLURM job submission with scalable job arrays
   - Experiment orchestration
3. Pre/post-processing layer (**MATLAB**)
   - Converts `.mat` forcing datasets into model-ready `.bin` inputs
   - Aggregates model output into `.mat` files
   - Computes derived diagnostics (e.g., biological pump metrics)

## Requirements

### Core model

- Fortran compiler (`ifort` or `gfortran`)
- Unix-like environment (Linux or macOS); HPC recommended for large or multi-column simulations
- [SLURM](https://slurm.schedmd.com/overview.html) workload manager (recommended for large ensemble or global simulations)

### Pre- and post-processing

- [MATLAB](https://mathworks.com/products/matlab.html) (R2021a or later)

### Optional MATLAB toolboxes

- Parallel Computing Toolbox (optional; used to accelerate bootstrap uncertainty calculations for BCP metrics. All routines run in serial mode if unavailable.)

### Third-party MATLAB utilities

The following packages from [MATLAB's File Exchange](https://mathworks.com/matlabcentral/fileexchange/) are required and should be added to the MATLAB path: `m_map`, `brewermap`, `subaxis`

## Quick Start

### Step 1 - Create a configuration

**1.1 Copy the namelist template**

```bash
cp ./config/namelist_template.txt ./config/namelist_MYRUN.txt
```
Replace `MYRUN` with your experiment name. 

**1.2 Choose the grid type**

Edit your new namelist and set:
```
CHOICE_GRID_DOMAIN=2 ! local column(s)
CHOICE_GRID_DOMAIN=1 ! global grid
```
If running **local columns** (specific locations), also create in `./config/`:
```
config_latitudes_MYRUN.txt 
config_longitudes_MYRUN.txt
```
These files define the coordinates of your local columns.

**1.3 Adjust configuration files**

- **Scientific parameters** &rarr; `namelist`
- **Workflow options** &rarr; `./code/bash/core/slams_workflow_config.sh` 
- **Root directory paths** &rarr; `./code/bash/core/slams_paths_config.sh`

### Step 2 - Prepare input and compile

```bash
./code/bash/create_model_input.sh MYRUN  
```

### Step 3 - Run the model

```bash
./code/bash/run_model_balanced.sh MYRUN 
```
(SLURM optional; controlled in `./code/bash/core/slams_workflow_config.sh`)

### Step 4 – Post-process

```bash
./code/bash/read_model_output.sh MYRUN 
```
Generates: `runsoutput.mat`, `bcpmetrics.mat` (optional).

### Step 5 – Visualise

To generate figures locally, ensure the following files are available:
```bash
./tests/MYRUN/modelruns/runsoutput.mat
./tests/MYRUN/modelruns/bcpmetrics.mat
./tests/MYRUN/modelinputdata/grid_run.mat
./tests/MYRUN/modelinputdata/waterColNumDepthLayers.txt
```
These files are produced during model execution and post-processing. Once available on your local machine, run the MATLAB plotting scripts in `./code/matlab/plotting/`.

## Reproducibility

The example configuration provided in this repository reproduces the LOCALTS6 experiment used in the manuscript. All required configuration files are included. External forcing datasets must be generated as described below.

## Repository Structure

 - `code/`
    - `fortran/`: core model
      - `src/`: source code and Makefile (*see "Scripts Overview"*)
      - `build/`: compiled files and executable
    - `bash/`: workflow orchestration scripts
        - `pipelines/`: user-facing run scripts
        - `core/`: shared workflow modules (paths, environment, etc.)
        - `maintenance/`: cleanup and recovery tools
        - `experiments/`: sensitivity studies tools
    - `matlab/`: pre- and post-processing tools
      - `entrypoints/`: Bash entry points
      - `modelinterface/`: Fortran mirror
      - `diagnostics/`: model-derived calculations
      - `io/`: file loading and processing
      - `plotting/`: figures
      - `datacompilation/`: observational datasets processing and plotting
      - `processformulations/`: biological/physical process formulations
      - `checks/`: formulation checks
- `config/`: experiment configuration (namelist, grid definitions)
- `modelresources/`: MATLAB utilities
    - `external/`: third-party functions (*see "Requirements" section*)
    - `internal/`: custom functions
- `data/`: input and derived datasets (*not provided*)
    - `raw/`: downloaded forcing datasets (*see below for details*) and observational compilation for model validation
    - `interim/`: intermediate forcing and pre-processed inputs
    - `processed/`: model-derived outputs and analysis products
- `tests/`: unit tests

## Forcing Data

Raw forcing data are not included due to licensing constraints.

Environmental forcing fields must be generated externally and placed in:
```bash
./data/raw/
```
The example configuration included in this repository was built using the following data products:
- Surface photosynthetic active radiation – NASA Aqua-MODIS sensor
- Net primary production (NPP) – BICEP project
- Chlorophyll *a* – ESA OC-CCI merged-sensor product
- Temperature, salinity, nutrients, dissolved oxygen – World Ocean Atlas 2023 (WOA23) 
- Aeolian dust flux – CMIP6/NCAR-CESM2
- Mesozooplankton biomass – CMIP6/PISCES
- Mixed layer depth – IFREMER
- Seawater density – calculated from WOA23 using GSW Oceanographic Toolbox
- Seawater dynamic viscosity – calculated from WOA23 using MIT seawater properties library routines
- Calcite saturation state – calculated using CO2SYS

Data download, formatting, and gap-filling tools are available in:
 - [ocean-data-lab](https://github.com/annarufas/ocean-data-lab)
 - [gap-filling-methods-ocean-data](https://github.com/annarufas/gap-filling-methods-ocean-data)

These tools convert raw NetCDF products into the `.mat` files required by SLAMS-2.0.

## Source Code

The SLAMS-2.0 model consists of the following Fortran 90 files:

| File name                     | Purpose                                                 |
|-------------------------------|---------------------------------------------------------|
| `main.F90`                    | Program entry point |
| `particlestructure.F90`       | Initialises the particle module |   
| `modelconstants.F90`          | Physical and chemical constants |
| `modelparameters.F90`         | Model parameters and configuration |
| `modelgrid.F90`               | Grid framework |
| `modelforcingdata.F90`        | Reads forcing inputs |
| `waterphysicsandlight.F90`    | Physical/light-related fields calculation |
| `initialisation.F90`          | Model initialisation |
| `injectsurfaceparticles.F90`  | Surface particle injection |
| `primaryproduction.F90`       | Primary production processes |
| `particledynamics.F90`        | Coagulation, breakup |
| `heterotrophicmetabolism.F90` | Zooplankton and bacterial metabolism |
| `mineraldissolution.F90`      | Abiotic mineral dissolution |
| `sink.F90`                    | Particle sinking, traps and imaging |
| `montecarlosampling.F90`      | Stochastic sampling routine |
| `calcparticleattributes.F90`  | Calculate particle attributes |
| `wrappers.F90`                | Module integration |
| `modelcounters.F90`           | Tracks diagnostic variables |
| `modeleulerianvariables.F90`  | Source-minus-sinks and auxiliary fields |
| `modelparticlecollection.F90` | Arrays to collect and classify particles for output |
| `modeloutput.F90`             | Output writer |
| `findfunctions.F90`           | Particle selection utilities |
| `sanitychecks.F90`            | Physical and numerical consistency checks |  
| `safemath.F90`                | Floating-point safety utilities |  
| `blockdefinitions.h`          | Defines block structures |
| `pstruct.h`                   | Particle attribute definitions |
| `timer.F90`                   | Simulation clock |
| `theseed.F90`                 | RNG control |

## Acknowledgments

SLAMS-2.0 originated as part of my PhD research at the University of Oxford (2016–2022) under the NERC large grant *COMICS* (Controls over Ocean Mesopelagic Interior Carbon Storage; NE/M020835/2). Additional support was provided by the University of Oxford’s COVID-19 Scholarship Extension Fund and Wolfson College COVID-19 Hardship Fund. Development continued with support from the European Space Agency (ESA) Living Planet Fellowship *SLAM DUNK* (Combining a Stochastic Lagrangian Model of Marine Particles with ESA’s Big Data to Understand the Effects of a Changing Ocean on the Planktonic Food Web; Contract No. 4000144464/24/I-DT-lr).

The original model concept was developed by [Dr. Tinna Jokulsdottir and Prof. David Archer](https://doi.org/10.5194/gmd-9-1455-2016), whose work laid the foundation for SLAMS-2.0.

If SLAMS-2.0 contributes to your research and you are interested in collaboration or supporting ongoing development, please feel free to contact me: anna.rufas@gmail.com.

**IMPORTANT**: Please do NOT post the SLAMS-2.0 code or any files downloaded from this repository on your own GitHub or other website. See LICENSE for licensing information.