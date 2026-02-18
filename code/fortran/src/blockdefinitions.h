! Use these to print information to the terminal screen
#undef BLOCK_PRINT_INFO_TO_THE_SCREEN
#undef BLOCK_PRINT_INFO_TO_THE_SCREEN_EXTENDED

! Use this to write information to output files
#define BLOCK_DIAGNOSTICS

! Use this to activate or deactivate sections of the code needed/not needed when coupling the model to the TMM
#undef BLOCK_TMM_SLAMS

! Use this to free up space in the particle array (not used by the TMM)
#define BLOCK_SQUEEZE_PARTICLE_ARRAY

! Use this to generate non-repeatable sequences of random numbers (i.e., reinitialise the 
! pseudo-random number generator). Leave it undef and, by default, the machine will give 
! rise to a repeatable sequence of random numbers (the desirable option).
! Read: https://stackoverflow.com/questions/23875589/why-are-my-random-numbers-always-the-same
#undef BLOCK_NEW_RANDOM_SEED

! Use this to save disk space when running with the optimiser
#define SAVE_DISK_SPACE

! Use these to debug the model by switching them off and then swithching them on gradually
#undef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
#undef BLOCK_FEW_TIME_STEPS
#define BLOCK_PHYTOPLANKTON
#define BLOCK_TEP
#define BLOCK_CLAY
#define BLOCK_COAGULATION
#define BLOCK_ZOOPLANKTON
#define BLOCK_MICROBIAL_METABOLISM
#define BLOCK_ABIOTIC_DISSOLUTION
#define BLOCK_ABIOTIC_FRAGMENTATION
#define BLOCK_FREE_LIVING_BACTERIA
#define BLOCK_SINK
#define BLOCK_PHOTOLYSIS

! seed
! add fragmentation (lead by fragmentation of non-sticky particles)
! add packing
! add photolysis
! add respiration
! add zoo
! add solub