#include "blockdefinitions.h"

module modelparameters

use modelconstants, only: PI, MESOZOO_DENSITY

implicit none
external :: findunit
public

! ----------------------------------------------------------------------------------------
! COMPILE-TIME STRUCTURAL CONSTANTS
!
! These parameters define memory layout, array sizes, and model structure.
! Changing any of these requires recompilation.
! ----------------------------------------------------------------------------------------

integer, parameter :: maxNumClusters            		= 5e5, &
					  nPfts                     		= 4,   & ! num. of phytoplankton functional types (diatoms, flgellated phyto, coccos, picophyto)				  
					  nMinerals                 		= 3,   & ! (opal, calcite, clay)
					  nSizeClasses              		= 16,  & ! as in EcoTaxa
					  nVeloClasses              		= 16,  &
					  nVolumeClasses            		= 27,  &
					  nMesoZooSizeClasses       		= 26,  &
					  nMainParticleTypes        		= 9,   & ! 1=faecal aggregate, 2=non-faecal aggregate, 3=living phyto, 4=dead phyto, 5=TEP, 6=clay, 7=dead zoo, 8=inorganic, single, non-faecal, 9=organic, single, non-faecal that is not dead phyto (fragment)
					  nTrackedParticles         		= 5,   & ! this is for purposes of checking model behaviour
					  maxNumSmsTerms            		= 50,  & ! sources-minus-sinks term 
					  maxNumTracers             		= 10,  & ! DOC, DIC, CO2, silicic acid, phosphate, nitrate, ALK
					  maxNumAuxTerms            		= 40,  & ! 
					  maxNumSedTrapDeployDepths 		= 58,  & ! for particle fluxes
					  maxNumImagingDeployDepths 		= 58,  & ! for particle spectra (make sure maxNumSedTrapDeployDepths = maxNumImagingDeployDepths)
					  nNewPhytoClusters         		= 19,  & ! 
					  nNewClayClusters          		= 1,   & ! 
					  nNewBactTepClusters       		= 1,   &
					  nNewPhytoTepClusters      		= 2,   & !
					  nNewFixedClusters         		= 1e4, & ! for fixed initial conditions (to reproduce Smoluchowski solution)
					  nNewZooDeadClustersPerDepthLayer 	= 1,   &
					  maxNumZooDeadClustersPerProfile  	= 9

! ----------------------------------------------------------------------------------------
! ENUMERATED INDEX CONSTANTS
!
! Symbolic integer identifiers used to index structured model arrays (e.g., Lagrangian
! particle composition, SMS terms, auxiliary diagnostics).
! ----------------------------------------------------------------------------------------

! Indexes in the Lagrangian particle variables
integer, parameter :: iOpal = 1, iCalcite = 2, iClay = 3
		  
! Row indexes in the SMSterm array (units of mol m-3 d-1)
integer, parameter :: iPrimProdOrgC     = 1,  & ! C taken up by phytoplankton
					  iPrimProdCaCO3    = 2,  & ! CO3= taken up by coccos
					  iPrimProdOpal     = 3,  & ! SiOH4 taken up by diatoms					  
					  iDepoClay	        = 4,  & ! aeolian dust deposition
					  iProdTepPhyto     = 5,  & ! TEC released by phytoplankton during photosynthesis
					  iProdTepMicrob    = 6,  & ! TEC released by microbes
					  iPhotoTepC        = 7,  & ! DOC released due to photolysis	  
					  iZooSolubOrgC     = 8,  & ! DOC released due to zooplankton messy feeding
					  iZooSolubTepC     = 9,  & ! dissolved TEC released due to zooplankton messy feeding
					  iZooSolubCaCO3    = 10, & ! CO3= released due to zooplankton messy feeding
					  iZooSolubOpal     = 11, & ! SiOH4 released due to zooplankton messy feeding
					  iZooSolubClay     = 12, & ! dissolved clay released due to zooplankton messy feeding
					  iZooIngestOrgC    = 13, & ! zooplankton ingestion of POC
					  iZooIngestTepC    = 14, & ! zooplankton ingestion of TEC
					  iZooRespOrgC      = 15, & ! C (CO2) released due to zooplankton respiration of POC
					  iZooRespTepC      = 16, & ! C (CO2) released due to zooplankton respiration of TEC
					  iZooEgestOrgC     = 17, & ! zooplankton egestion of POC
					  iZooEgestTepC     = 18, & ! zooplankton egestion of TEC
					  iZooDissolCaCO3   = 19, & ! CO3= release due to biotic dissolution of CaCO3
					  iZooExcretOrgC    = 20, & ! DOC release due to zooplankton excretion
					  iZooExcretTepC    = 21, & ! dissolved TEC release due to zooplankton excretion
					  iZooDeathOrgC  	= 22, & ! POC released as zooplankton dead bodies
					  iZooDeathCaCO3 	= 23, & ! CaCO3 released as zooplankton dead bodies
					  iZooDeathOpal  	= 24, & ! opal released as zooplankton dead bodies
					  iMicrobRespOrgC   = 25, & ! C (CO2) released due to microbial respiration of POC
					  iMicrobRespTepC   = 26, & ! C (CO2) released due to microbial respiration of TEC
					  iMicrobSolubOrgC  = 27, & ! DOC released due to microbial solubilisation
					  iMicrobSolubTepC	= 28, & ! dissolved TEC released due to microbial solubilisation
					  iMicrobSolubCaCO3 = 29, & ! CO3= released due to abiotic solubilisation
					  iMicrobSolubOpal  = 30, & ! SiOH4 released due to abiotic solubilisation
					  iMicrobSolubClay  = 31, & ! dissolved clay released due to abiotic solubilisation					  
					  iDissolCaCO3      = 32, & ! CO3= released due to abiotic dissolution of CaCO3	
					  iDissolOpal       = 33    ! SiOH4 released due to abiotic dissolution of opal
					  			  
! Row indexes in the auxTerm array
integer, parameter :: iZooBiomass             		= 1,  & ! mol C m-3
					  iZooNumber              		= 2,  & ! ind. m-3
					  iZooDeadBiomass         		= 3,  & ! mol C m-3
					  iDiatBiomass            		= 4,  & ! mol C m-3
					  iFlagelBiomass                = 5,  & ! mol C m-3
					  iCoccoBiomass          		= 6,  & ! mol C m-3
					  iPicoBiomass            		= 7,  & ! mol C m-3
					  iFreshDiatCellQuota     		= 8,  & ! mol C cell-1
					  iFreshFlagelCellQuota     	= 9,  & ! mol C cell-1
					  iFreshCoccoCellQuota    		= 10, & ! mol C cell-1
					  iFreshPicoCellQuota     		= 11, & ! mol C cell-1				  
					  iAvgProbDiat	                = 12, &
					  iAvgProbFlagel	            = 13, &
					  iAvgProbCocco	                = 14, &
					  iAvgProbPico	                = 15, &	
					  iProfProbDiat                 = 16, &
					  iProfProbFlagel               = 17, &
					  iProfProbCocco                = 18, &
					  iProfProbPico                 = 19, &
					  iZooSpecRespRate        		= 20, & ! s-1
					  iMicrobSpecRespRate     		= 21, & ! s-1
					  iNightDvmUpperBound     		= 22, & ! m
					  iNightDvmLowerBound     		= 23, & ! m
					  iDayDvmUpperBound       		= 24, & ! m
					  iDayDvmLowerBound	      		= 25, & ! m
					  iEuphoticDepth          		= 26, & ! m
					  iCollisionKernel        		= 27, & ! m3 s-1
					  iBrownianKernel	      		= 28, & ! m3 s-1
					  iShearKernel	          		= 29, & ! m3 s-1
					  iSettlingKernel	      		= 30, & ! m3 s-1
					  iNumClusterPairsEvalCoagu 	= 31, & ! num.
					  iNumPxCmore		        	= 32, &
					  iNumPxCless		        	= 33, &
					  iCoagulationSuccess     		= 34, &
					  iCoagulationProbability 		= 35, &
					  iNumParticlesEvalEncounter   	= 36, &
					  iNumParticlesEncountered	   	= 37, &	
					  iNumParticlesFragmentedZoo    = 38, &
					  iNumParticlesIngestedZoo     	= 39, &
					  iNumParticlesOmittedZoo    	= 40

! ----------------------------------------------------------------------------------------
! GLOBAL MODEL CONSTANTS (NON-STRUCTURAL)
!
! Non-structural compile-time constants defining global modelling assumptions and
! numerical resolution settings inherited from SLAMS-1.0.
!
! These do NOT affect memory layout but require recompilation if changed.
! ----------------------------------------------------------------------------------------

real*8, parameter :: maxPrimParticlesPerAggregateForCollision = 1d9,    & 
					 decreaseCollisionTimeframeFactor         = 0.10d0, &
					 timeStep                                 = 28800d0 ! s (num of seconds in 1/3 of a day, =24d0*3600d0/3d0) 28800

! ----------------------------------------------------------------------------------------
! PERSISTENT MODEL STATE ARRAYS
!
! Static arrays whose dimensions depend on compile-time structural constants.
! ----------------------------------------------------------------------------------------

real*8, dimension(nSizeClasses) :: particleSizeClasses
real*8, dimension(nVeloClasses) :: particleVeloClasses
real*8, dimension(nVolumeClasses) :: particleVolumeClasses
real*8, dimension(nMesoZooSizeClasses+1) :: zooSizeClassLimits, zooWetWeightQuotaLimits
real*8, dimension(nMesoZooSizeClasses) :: zooSizeClasses, zooRadius, zooVolume, zooWetWeightQuota, &
	zooCarbonWeightQuota, zooSwimmingSpeed, zooProsomeLength, zooWetWeightQuotaBinWidth
real*8 :: sumZooWetWeightProducts
real*8, dimension(maxNumSedTrapDeployDepths) :: possibleSedTrapDeployDepths
real*8, dimension(maxNumImagingDeployDepths) :: possibleImagingDeployDepths
									  
! ----------------------------------------------------------------------------------------
! RUNTIME-CONFIGURABLE SCIENTIFIC PARAMETERS (NAMELIST)
!
! These parameters are read from namelist.input and can be modified between runs without 
! recompiling the model.
!
! They define the scientific experiment configuration.
! ----------------------------------------------------------------------------------------

! Reference biological carbon pump parameters
real*8 :: fractal_dimension_agg_max, fractal_dimension_agg_min, porosity_faecalpell_max, &
	porosity_faecalpell_min, stickiness_phyto_cell_initial
real*8 :: C_quota_diat_max, C_quota_diat_min, C_quota_flagel_max, C_quota_flagel_min, &
	C_quota_cocco_max, C_quota_cocco_min, C_quota_pico_max, C_quota_pico_min, Si2C_diat, &
	Calc2C_cocco_max, k_omega, phyto_life_span_factor, phyto_max_allowed_life_days
real*8 :: C_quota_TEP_max, C_quota_TEP_min, phyto_exudation_frac, k_TEP
real*8 :: zoo_distrib_slope, detection_radius_factor_mesozoo, agg_to_zoo_size_ratio, &
	frac_OM_zoo_ingestion_surf, frac_OM_zoo_ingestion_deep, carbon_density_threshold_zoo_ingestion, &
	zoo_absorption_eff_carbon, zoo_net_growth_eff, gut_passage_time_mesozoo, mort_rate_mesozoo 
real*8 :: q10_mesozoo, q10_microb, resp_rate_poc_max_0deg_mesozoo, resp_rate_poc_max_0deg_microb, &
	resp_rate_tepc_max_0deg_microb, resp_rate_poc_max_mesozoo, resp_rate_poc_max_microb, &
	resp_rate_tepc_max_microb, k_O2_resp, solub_rate_poc, solub_rate_tepc, photodegradation_rate_tepc
real*8 :: dissol_rate_calc, dissol_rate_calc_zoo_gut, dissol_rate_opal_0deg, q10_bSi

! Fixed parameters
real*8 :: sea_surface_microlayer_depth, breaking_reynolds_threshold, min_par_for_photosynthesis, &
	order_react_dissol_calc, detection_limit_poc, detection_limit_calc, detection_limit_opal, &
	detection_limit_clay, operational_size_poc_min, operational_size_doc_min, growth_rate_max_0deg_diat, &
	growth_rate_max_0deg_flagel, growth_rate_max_0deg_cocco, growth_rate_max_0deg_pico, &
	growth_rate_max_phyto, alpha_chl_spec_diat, alpha_chl_spec_flagel, alpha_chl_spec_cocco, &
	alpha_chl_spec_pico, q10_diat, q10_flagel, q10_cocco, q10_pico, k_NO3_diat, k_NO3_flagel, &
	k_NO3_cocco, k_NO3_pico, k_PO4_diat, k_PO4_flagel, k_PO4_cocco, k_PO4_pico, k_Si_diat, &
	Chl2C_ratio, C2N_ratio, C2O_ratio, C2P_ratio, C_frac_in_OM, C_frac_in_TEP, clay_quota
	
! Configuration parameters 
integer :: nDepthLayers, nYears, nDaysPassedToShiftClusters, nDaysPassedToSeedClays, &
	nDaysPassedToDestroyLargeAggs
integer :: choicePhytoSeedingDepthDistrib, choicePhytoIrradLimFunc, choiceIsPhytoGrowthLim, &
	choicePhytoCellKillingScheme, choiceFractalDimensionScheme, choiceStickinessFunc, &
	choiceMineralEffOnFaecalPellPorosity, choiceCollisionSamplingScheme, choiceCollisionKernels, &
	choiceIncludeHydrodynamicForces, choiceCriteriaToBreak, shrinkAfterMineralDissolutionScheme, &
	shrinkAfterTepDegradationScheme, shrinkAfterMicrobialMetabolismScheme, choiceZooBehaviourKernels, &
	choiceZooIngestionCriteria
logical :: choiceIsSinkingAffectedByTurbulence, choiceIsMineralProtectAgainstMicrobResp

! ----------------------------------------------------------------------------------------
! DERIVED RUNTIME STATE VARIABLES
!
! Computed from the full runtime configuration (scientific parameters, numerical 
! resolution, and global constraints) to ensure internal consistency.
! ----------------------------------------------------------------------------------------

integer :: nTimeStepsDay,nTimeStepsYear, maxNumTimeSteps, nTimeStepsPassedToShiftClusters, &
	nTimeStepsPassedToSeedClay, nTimeStepsPassedToDestroyLargeAggs, maxNumClustersPerProfile
real*8, dimension(nPfts) :: C_quota_pft_max, C_quota_pft_min, growth_rate_max_0deg_pft, &
	k_NO3_pft, k_PO4_pft, alpha_chl_spec_pft, q10_phyto

! ----------------------------------------------------------------------------------------
! RUNTIME INTERFACE DEFAULTS (NO NAMELIST OVERRIDE)
!
! Deterministic binary filenames defining the interface between:
!   - MATLAB preprocessing and post-processing
!   - The Fortran model core
!
! These are intentionally NOT user-configurable.
! If new forcing variables or diagnostics are introduced, this section and the 
! corresponding MATLAB routines must be updated consistently.
! ----------------------------------------------------------------------------------------

character(len=*), parameter :: &
    filenamePAR0      				= 'in_par0.bin', &
    filenameNPP       				= 'in_npp.bin', &
    filenameChla      				= 'in_chla.bin', &
    filenameRho       				= 'in_rho.bin', &
    filenameOmegaCalc 				= 'in_omegacalcite.bin', &
    filenameNO3       				= 'in_nitrate.bin', &
    filenameSiOH4     				= 'in_silicicacid.bin', &
    filenamePO4       				= 'in_phosphate.bin', &
    filenameO2        				= 'in_oxygen.bin', &
    filenameTempC     				= 'in_temperature.bin', &
    filenameAeolClay  				= 'in_dust.bin', &
    filenameMLD       				= 'in_mld.bin', &
    filenameMesoZoo   				= 'in_mesozooplankton.bin', &
    filenameDynVisco  				= 'in_dynamicviscosity.bin', &
    filenameZub       				= 'in_depthupperbounds.bin', &
    filenameZlb       				= 'in_depthlowerbounds.bin'

character(len=*), parameter :: &
	filenameFlux                    = 'out_flux.bin', &
	filenameFluxSf                  = 'out_flux_seafloor.bin', &
 	filenameAvgAttSizeClass         = 'out_avgparticleatt_sizeclass.bin', &
	filenameAvgAttVeloClass			= 'out_avgparticleatt_veloclass.bin', &
	filenameAvgAttVolumeClass		= 'out_avgparticleatt_volumeclass.bin', &
	filenameAvgAttSizeClassSf		= 'out_avgparticleatt_sizeclass_seafloor.bin', &
	filenameAvgAttVeloClassSf		= 'out_avgparticleatt_veloclass_seafloor.bin', &
	filenameInstAvgAttSizeClass		= 'out_inst_avgparticleatt_sizeclass.bin', &
	filenameInstAvgAttVeloClass		= 'out_inst_avgparticleatt_veloclass.bin', &
	filenameInstAvgAttSizeClassSf	= 'out_inst_avgparticleatt_sizeclass_seafloor.bin', &
	filenameInstAvgAttVeloClassSf	= 'out_inst_avgparticleatt_veloclass_seafloor.bin', &
	filenameAvgAttMainType			= 'out_avgparticleatt_maintype.bin', &
	filenameAvgAttMainTypeSf		= 'out_avgparticleatt_maintype_seafloor.bin', &
	filenameInstAvgAttMainType		= 'out_inst_avgparticleatt_maintype.bin', &
	filenameInstAvgAttMainTypeSf	= 'out_inst_avgparticleatt_maintype_seafloor.bin', &
	filenameLossTerms				= 'out_lossterms.bin', &
	filenameSms						= 'out_sms.bin', &
	filenameSmsIntegrated			= 'out_sms_integrated.bin', &
	filenameAuxTerms				= 'out_aux.bin', &
	filenameStatsNumClusters		= 'out_stats_nClusters.bin', &
	filenameStatsNumParticles		= 'out_stats_nParticles.bin', &
	filenameControl					= 'out_control.bin', &
	filenameClustersInteger			= 'out_cluster_i.bin', &
	filenameClustersReal            = 'out_cluster_r.bin'
	
contains

! ========================================================================================

subroutine EstablishParticleSizeAndVeloCategories( )

	integer :: iSc, iVc, iVolClass, i
	character(len=30) :: fmt

	particleSizeClasses(:) = 0d0
	particleVeloClasses(:) = 0d0	
	
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS

	! As in 'Simple' (G. Jackson) - those are the upper boundaries
	particleVolumeClasses(:) = 0d0
	do iVolClass = 1, nVolumeClasses
		particleVolumeClasses(iVolClass) = 1d12*((PI/6d0)*(10d-4)**3) * 2d0**iVolClass ! cm3 -> um3
	end do
	
#endif

	! As in EcoTaxa
	! 16 diameter limits (um) : 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1020, 2050, 4100, 8190, 16400, 32768
	! The first diameter bin is 0.2–1 um (including 1 um), and the last bin is >16400 um

	do iSc = 1, nSizeClasses
		if (iSc == 1) then
			particleSizeClasses(iSc) = 1d0
		else
			particleSizeClasses(iSc) = 2d0*particleSizeClasses(iSc-1) ! geometric sequence
		end if
	end do

	! Similarly used in Trull et al. (2008)
	! 16 velocity limits (m d-1) : 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102, 205, 410, 819, 1640, >1640
	! The first velocity bin is 0–0.1 m d-1, and the last bin is 1640-3277 m d-1
	
	do iVc = 1, nVeloClasses
		if (iVc == 1) then
			particleVeloClasses(iVc) = 0.1d0
		else
			particleVeloClasses(iVc) = 2d0*particleVeloClasses(iVc-1) ! geometric sequence
		end if
	end do
	
	fmt = '(5F12.1)'	
	write(*,*)
	write(*,*) 'Particle size (diameter) classes (um):'
	write(*,fmt) (particleSizeClasses(i), i=1,nSizeClasses)
	write(*,*)
	write(*,*) 'Particle velocity classes (m d-1):'
	write(*,fmt) (particleVeloClasses(i), i=1,nVeloClasses)
	write(*,*)
	
end subroutine EstablishParticleSizeAndVeloCategories

! ========================================================================================

subroutine EstablishMesoZooSizeCategoriesAndSpectraVariables( )

	integer :: iZc
	real*8, dimension(nMesoZooSizeClasses+1) :: zooRadiusLimits, zooVolumeLimits, &
		zooCarbonWeightQuotaLimits

	zooSizeClassLimits(:) = 0d0
	zooSizeClasses(:) = 0d0
	
	zooRadius(:) = 0d0
	zooRadiusLimits(:) = 0d0
	
	zooVolume(:) = 0d0
	zooVolumeLimits(:) = 0d0
	
	zooCarbonWeightQuota(:) = 0d0
	zooCarbonWeightQuotaLimits(:) = 0d0
	
	zooWetWeightQuota(:) = 0d0
	zooWetWeightQuotaLimits(:) = 0d0
	zooWetWeightQuotaBinWidth(:) = 0d0
	
	zooSwimmingSpeed(:) = 0d0
	zooProsomeLength(:) = 0d0
	
	! This size sequence is similar to Kwong & Pakhomov (2021), Suppl. Table S1. It is based on 
	! equivalent spherical diameter (ESD)
	do iZc = 1, nMesoZooSizeClasses+1
    	if (iZc == 1) then
        	zooSizeClassLimits(iZc) = 185d0 ! um
    	else
        	zooSizeClassLimits(iZc) = NINT(1.17d0*zooSizeClassLimits(iZc-1),8) ! um
    	end if
	end do

	! Geometric mean (middle point)
	do iZc = 1, nMesoZooSizeClasses
		zooSizeClasses(iZc) = (zooSizeClassLimits(iZc)*zooSizeClassLimits(iZc+1))**0.5d0 ! um
	end do

	! Swimming speed, based on Kiorboe 2011, Fig. 3
	zooSwimmingSpeed(:) = (10d0**(0.39d0 + (0.79d0*LOG10((zooSizeClasses(:)*1d-4)))))*1e-2 ! cm s-1 --> m s-1 (zooSizeClasses in cm)

	! Prosome length (ESD = PL x aspect_ratio^(2/3)), where aspect ratio is mentioned in Kiorboe (2011)
	zooProsomeLength(:) = zooSizeClasses(:)/0.38d0**(2d0/3d0)

	! Radius is half of ESD
	zooRadius(:) = zooSizeClasses(:)/2d0 ! um

	! Volume, based on Brun et al. 2019
	zooVolume(:) = (4d0/3d0)*PI*(zooRadius(:)*1d-4)**3 ! cm3 (radius in cm)

	! Wet weight quota, based on an assumed copepod density (as used in Brun et al. 2019, suggested by Kiorboe 2013)
	zooWetWeightQuota(:) = (zooVolume(:)*MESOZOO_DENSITY)*1d3 ! g WW --> mg WW
	! calculating it from volume or as the geometric mean from the weight wet quota limits yields the same result
	
	! Carbon weight quota, based on Kiorboe 2013, Table 2
	zooCarbonWeightQuota(:) = 10d0**(-0.93d0+(0.95d0*LOG10(zooWetWeightQuota(:)))) ! mg C

	! Class limits
	zooRadiusLimits(:) = zooSizeClassLimits(:)/2d0 ! um
	zooVolumeLimits(:) = (4d0/3d0)*PI*(zooRadiusLimits(:)*1d-4)**3	! cm3
	zooWetWeightQuotaLimits(:) = (zooVolumeLimits(:)*MESOZOO_DENSITY)*1d3 ! g WW --> mg WW
	zooCarbonWeightQuotaLimits(:) = 10d0**(-0.93d0+(0.95d0*LOG10(zooWetWeightQuotaLimits(:)))) ! mg C

	! Other spectrum parameters
	do iZc = 1, nMesoZooSizeClasses
    	zooWetWeightQuotaBinWidth(iZc) = zooWetWeightQuotaLimits(iZc+1)-zooWetWeightQuotaLimits(iZc)
	end do
    sumZooWetWeightProducts = SUM((zooWetWeightQuota(:)**(-zoo_distrib_slope))*zooWetWeightQuotaBinWidth(:)) ! Samar calls this 'g'
    
    !print *, 'zoo radius', zooRadius(:)
    !print *, 'zoo swimming speed', zooSwimmingSpeed(:)
    !print *, 'sumZooWetWeightProducts', sumZooWetWeightProducts
    !print *, 'zooWetWeightQuotaLimits', zooWetWeightQuotaLimits(:)

end subroutine EstablishMesoZooSizeCategoriesAndSpectraVariables

! ========================================================================================

subroutine GetDepthsForParticleDataCollection( )

	! The first deployment depth and imaging depth have to be shallow enough so that we can 
	! register the properties of shallow water columns (at least 20 m shallow)	

	real*8, dimension(maxNumSedTrapDeployDepths) :: depths ! remember to tweak maxNumSedTrapDeployDepths at the beginning of this script to catch all the depth intervals
	real*8 :: currentDepth
    integer :: idx
    
    possibleSedTrapDeployDepths = 0d0
    possibleImagingDeployDepths = 0d0
    depths = 0d0

    idx = 1
    currentDepth = 10d0

    ! 10 m intervals from 10 to 200 m
	do while (currentDepth <= 200d0 .and. idx <= maxNumSedTrapDeployDepths)
		depths(idx) = currentDepth
		currentDepth = currentDepth + 10d0
		idx = idx + 1
	end do

	! 50 m intervals from 250 to 1000 m
	currentDepth = 250d0
	do while (currentDepth <= 1000d0 .and. idx <= maxNumSedTrapDeployDepths)
		depths(idx) = currentDepth
		currentDepth = currentDepth + 50d0
		idx = idx + 1
	end do

	! 100 m intervals from 1100 to 2000 m
	currentDepth = 1100d0
	do while (currentDepth <= 2000d0 .and. idx <= maxNumSedTrapDeployDepths)
		depths(idx) = currentDepth
		currentDepth = currentDepth + 100d0
		idx = idx + 1
	end do

	! 250 m intervals from 2250 to 5000 m
	currentDepth = 2250d0
	do while (currentDepth <= 5000d0 .and. idx <= maxNumSedTrapDeployDepths)
		depths(idx) = currentDepth
		currentDepth = currentDepth + 250d0
		idx = idx + 1
	end do

	! Assign to public arrays
	possibleSedTrapDeployDepths = depths
	possibleImagingDeployDepths = depths

end subroutine GetDepthsForParticleDataCollection

! ========================================================================================

subroutine InitialiseRuntimeParameters()

	integer :: ioUnit, ioStatusRef, ioStatusFix, ioStatusConf
	
	! ------------------------------------------------------------------------------------
	! NAMELIST DECLARATIONS
	!
	! Define which runtime-configurable variables may be overridden by the user.
	! ------------------------------------------------------------------------------------
	
	namelist /ReferenceParameters/ &
		fractal_dimension_agg_max, fractal_dimension_agg_min, porosity_faecalpell_max, &
		porosity_faecalpell_min, stickiness_phyto_cell_initial, C_quota_diat_max, &
		C_quota_diat_min, C_quota_flagel_max, C_quota_flagel_min, C_quota_cocco_max, &
		C_quota_cocco_min, C_quota_pico_max, C_quota_pico_min, Si2C_diat, Calc2C_cocco_max, &
		k_omega, phyto_life_span_factor, phyto_max_allowed_life_days, C_quota_TEP_max, &
		C_quota_TEP_min, phyto_exudation_frac, k_TEP, zoo_distrib_slope, detection_radius_factor_mesozoo, &
		agg_to_zoo_size_ratio, frac_OM_zoo_ingestion_surf, frac_OM_zoo_ingestion_deep, &
		carbon_density_threshold_zoo_ingestion, zoo_absorption_eff_carbon, zoo_net_growth_eff, &
		gut_passage_time_mesozoo, mort_rate_mesozoo, q10_mesozoo, q10_microb, resp_rate_poc_max_0deg_mesozoo, &
		resp_rate_poc_max_0deg_microb, resp_rate_tepc_max_0deg_microb, resp_rate_poc_max_mesozoo, &
		resp_rate_poc_max_microb, resp_rate_tepc_max_microb, k_O2_resp, solub_rate_poc, &
		solub_rate_tepc, photodegradation_rate_tepc, dissol_rate_calc, dissol_rate_calc_zoo_gut, &
		dissol_rate_opal_0deg, q10_bSi
		
	namelist /FixedParameters/ &
		sea_surface_microlayer_depth, breaking_reynolds_threshold, min_par_for_photosynthesis, &
		order_react_dissol_calc, detection_limit_poc, detection_limit_calc, detection_limit_opal, &
		detection_limit_clay, operational_size_poc_min, operational_size_doc_min, &
		growth_rate_max_0deg_diat, growth_rate_max_0deg_flagel, growth_rate_max_0deg_cocco, &
		growth_rate_max_0deg_pico, growth_rate_max_phyto, alpha_chl_spec_diat, alpha_chl_spec_flagel, &
		alpha_chl_spec_cocco, alpha_chl_spec_pico, q10_diat, q10_flagel, q10_cocco, q10_pico, &
		k_NO3_diat, k_NO3_flagel, k_NO3_cocco, k_NO3_pico, k_PO4_diat, k_PO4_flagel, k_PO4_cocco, &
		k_PO4_pico, k_Si_diat, Chl2C_ratio, C2N_ratio, C2O_ratio, C2P_ratio, C_frac_in_OM, C_frac_in_TEP, &
		clay_quota
			
	namelist /ConfigurationParameters/ &
		nDepthLayers, nYears, nDaysPassedToShiftClusters, nDaysPassedToSeedClays, &
		nDaysPassedToDestroyLargeAggs, choicePhytoSeedingDepthDistrib, choicePhytoIrradLimFunc, &
		choiceIsPhytoGrowthLim, choicePhytoCellKillingScheme, choiceFractalDimensionScheme, &
		choiceStickinessFunc, choiceMineralEffOnFaecalPellPorosity, choiceCollisionSamplingScheme, &
		choiceCollisionKernels, choiceIncludeHydrodynamicForces, choiceCriteriaToBreak, &
		choiceIsSinkingAffectedByTurbulence, shrinkAfterMineralDissolutionScheme, &
		shrinkAfterTepDegradationScheme, shrinkAfterMicrobialMetabolismScheme, &
		choiceIsMineralProtectAgainstMicrobResp, choiceZooBehaviourKernels, choiceZooIngestionCriteria
		
	! ------------------------------------------------------------------------------------
	! DEFAULT RUNTIME VALUES
    !
    ! These defaults are applied before reading namelist.input.
    ! If the namelist provides values, they overwrite these.
	! ------------------------------------------------------------------------------------

	! ------ REFERENCE parameters defaults ------
	
	! Particle parameters
	fractal_dimension_agg_max 		= 1.9d0 !	look at Kilps (1993)'s thesis (for marine aggregates (>O.5mm) D3 has been found to be <2 (Logan & Wilkinson 1990))
	fractal_dimension_agg_min 		= 1.3d0
	porosity_faecalpell_max   		= 0.50d0
	porosity_faecalpell_min   		= 0.20d0
	stickiness_phyto_cell_initial 	= 0.50d0
	
	! Phytoplankton particle parameters
	C_quota_diat_max   			= 5000d-12 ! mol C cell-1
	C_quota_diat_min   			= 40d-12 ! mol C cell-1
	C_quota_flagel_max 			= 40000d-12 ! mol C cell-1
	C_quota_flagel_min 			= 9000d-12 ! mol C cell-1
	C_quota_cocco_max  			= 20d-12 ! mol C cell-1
	C_quota_cocco_min  			= 5.0d-12 ! mol C cell-1
	C_quota_pico_max   			= 1.0d-12 ! mol C cell-1
	C_quota_pico_min   			= 0.080d-12! mol C cell-1	
	Si2C_diat          			= 0.20d0 ! based on ~1 mol Si:1 mol N
	Calc2C_cocco_max   			= 1.04d0 ! from Gangsto et al. (2011)
	k_omega      				= 0.40d0 
	phyto_life_span_factor		= 2d0 ! typical default (cells live up to ~2 doublinf times)
	phyto_max_allowed_life_days = 30d0

	! TEP particle parameters	
	C_quota_TEP_max      = 6000d-12 ! mol C TEP-1
	C_quota_TEP_min    	 = 3.0d-12	! mol C TEP-1	
	phyto_exudation_frac = 0.20d0
	k_TEP                = 0.10d0 ! TEP fraction at stickiness=0.5
	
	! Zooplankton parameters	 			
	zoo_distrib_slope                		= 0.80d0 ! Kwong & Pakhomov (2021)
	detection_radius_factor_mesozoo  		= 3d0
	agg_to_zoo_size_ratio            		= 2d0
	frac_OM_zoo_ingestion_surf       		= 0.30d0
	frac_OM_zoo_ingestion_deep       		= 0.30d0
	carbon_density_threshold_zoo_ingestion 	= 1d-6 ! g C cm-3
	zoo_absorption_eff_carbon        		= 0.69d0
	zoo_net_growth_eff               		= 0.75d0	
	gut_passage_time_mesozoo         		= 4000d0 ! s (standard value is 30 min = 1,800 s, Irigoien (1998))
	mort_rate_mesozoo                		= 0.020d0 ! d-1

	! Q10 factors					
	q10_mesozoo  = 2.0d0 ! 1.9, Ikeda et al. 2014
	q10_microb   = 1.65d0 ! Q10-factor within the range of 2–3 in large areas of the ocean [Pomeroy and Wiebe, 2001],
	
	! Respiration rates, carbon-specific, d-1 (fraction of body carbon respired per day) 
	resp_rate_poc_max_0deg_mesozoo  = 0.050d0 ! Hernandez-Leon & Ikeda, 2005
	resp_rate_poc_max_0deg_microb   = 0.050d0 ! 0.05 (Bendtsen et al. 2018)
	resp_rate_tepc_max_0deg_microb  = 0.09d0
	resp_rate_poc_max_mesozoo       = 0.25d0	
	resp_rate_poc_max_microb        = 0.20d0 ! maximum physiological respiration rate of bacteria (Cavan et al. 2017)
	resp_rate_tepc_max_microb       = 0.50d0 ! Mari et al. (2017)

	! Other degradation rate parameters
	k_O2_resp                  = 0.0045 ! mL L-1 (Morris & Schmidt, 2013, García-Robledo et al., 2016) (Gong et al. 2015 suggest 20-200 nmol L-1, which (our current value equals 200 nmol L-1))
	solub_rate_poc             = 0.010d0 ! d-1
	solub_rate_tepc            = 0.010d0 ! d-1
	photodegradation_rate_tepc = 0.30d0 ! d-1, Mari et al. (2017)
	
	! Abiotic mineral dissolution parameters
	dissol_rate_calc 		 = 5.0d0 ! d-1, kappa, Jansen et al. 2002 (range 1-7 d-1) (more stable than aragonite, so lower dissolution rate constant)
	dissol_rate_calc_zoo_gut = 20d0 ! d-1, Jansen & Wolf Gladrow (2001)
	!dissol_rate_opal = 0.05d0 ! d-1 (0.05 approx. d-1 from VanCappellen et al. 2002, World average surface waters; 0.008 to 0.18 from Bidle & Azam 1999, in this case when colonized by bacteria)
	dissol_rate_opal_0deg 	 = 1d-3
	q10_bSi 				 = 2.3d0 ! Kamatani 1982

	! ------ FIXED parameters defaults ------

	sea_surface_microlayer_depth = 1d-4 ! m (Mari et al. 2017)
	breaking_reynolds_threshold  = 20d0 ! 20-30 (Alldredge & Gotschalk 1988; Kiorboe et al. 2001)
	min_par_for_photosynthesis   = 1d-3 ! photosynthesis is possible where PARz > 0.1% of PAR0 (Buesseler et al. 2020)
	order_react_dissol_calc = 2d0 ! (Clancy recommends value of 2) (from 1 - 4.5, Sarmiento & Gruber)

	! Below what absolute amount does this quantity becomes meaningless for the model?
	! From flow cytometry and NanoSIMS (Secondary Ion Mass Spectrometry) analysis, detection 
	! limit is at ~1 fg C per particle = 8.3 x 10-17 mol C per particle.
	! From SEM–EDS, detection limits for minerals are: 
	!	5 fg CaCO3 = 5.0 x 10-17 mol per particle
	!	10 fg opal = 1.5 x 10-16 mol per particle
	!	50 fg clay = 1.3 x 10-16 mol per particle
	detection_limit_poc  = 8.3d-17 ! mol C 
	detection_limit_calc = 5.0d-17 ! mol CaCO3 
	detection_limit_opal = 1.5d-16 ! mol opal 
	detection_limit_clay = 1.3d-16 ! mol clay 

	operational_size_poc_min = 0.20d0 ! um (Kharbush et al. 2020, diameter)
	operational_size_doc_min = 1d-3 ! um (=1 nm)
	
	! Maximum specific growth rate at 0ºC (ranges), day-1, own compilation
	! A specific growth rate of 1 day-1 equals (1/0.70) 1.4 doublings day-1 
	! There is an important relationship between the percent growth rate and its 
	! doubling time known as “the rule of 70”: to estimate the doubling time for a 
	! steadily growing quantity, simply divide the number 70 by the percentage growth rate. 
	! Note: biomass-specific growth rates are expressed in mg C m-3 d-1
	growth_rate_max_0deg_diat   = 0.70d0
	growth_rate_max_0deg_flagel = 0.53d0 ! smaller than cocco
	growth_rate_max_0deg_cocco  = 0.38d0
	growth_rate_max_0deg_pico   = 0.35d0
	growth_rate_max_phyto       = 1.20d0 ! Flynn et al. 2017

	! Chla-specific initial slope of photosynthesis-irradiance curve (i.e. ratio PI at 
	!very low light), g C (g Chla)-1 m2 (umol photons)-1), from Geider et al. 1997		
	alpha_chl_spec_diat   = 3.7d-5
	alpha_chl_spec_flagel = 0.48d-5
	alpha_chl_spec_cocco  = 0.30d-5
	alpha_chl_spec_pico   = 2.5d-5

	! Q10 factors	
	q10_diat     = 2.2d0
	q10_flagel   = 2.08d0
	q10_cocco    = 1.88d0
	q10_pico     = 1.88d0
	
	! Phytoplankton nutrient half-saturation constants, mmol m-3 (=umol L-1), own compilation			
	k_NO3_diat   = 0.65d0
	k_NO3_flagel = 0.50d0
	k_NO3_cocco  = 0.40d0
	k_NO3_pico   = 0.23d0
	k_PO4_diat   = 0.05d0 ! based on Krumhardt et al. 2019 and Tyrrell & Taylor 1996
	k_PO4_flagel = 0.08d0
	k_PO4_cocco  = 0.005d0
	k_PO4_pico   = 0.004d0
	k_Si_diat    = 2.0d0 ! Balch et al. 2016

	! Other material parameters		
	Chl2C_ratio              = 2.5d-2 ! =40 g C (g Chl)-1 (max is 5.0d-2 ! =20 g C (g Chl)-1 (PISCES & TOPAZ), min is 5.0d-3 ! =200 g C (g Chl)-1) (Arteaga et al. 2016; Cloern 1995; Geider 1987)
	C2N_ratio                = 106d0/16d0 ! mol C (mol N)-1
	C2O_ratio                = 106d0/138d0 ! mol C (mol O2)-1
	C2P_ratio                = 106d0/1d0 ! mol C (mol P)-1
	C_frac_in_OM             = 0.54d0
	C_frac_in_TEP            = 0.95d0 ! fraction of carbon in tep
	clay_quota               = 1.0d-13
	
	! ------ CONFIGURATION parameter defaults ------

	nDepthLayers                  		  = 408
	nYears                        		  = 10 ! run duration, yr	
	nDaysPassedToShiftClusters            = 5
	nDaysPassedToSeedClays                = 10	
	nDaysPassedToDestroyLargeAggs         = 3
	
	choicePhytoSeedingDepthDistrib        	= 2 ! 1=SLAMS1.0 (fixed depth distribution), 2=SLAMS2.0 (depth distribution that is environment dependent)
	choicePhytoIrradLimFunc               	= 3 ! 1=Michaelis-Menten, 2=Exponential, 3=Smith, 4=Steele, 5=Tangent, 6=Plhotoinh.
	choiceIsPhytoGrowthLim                	= 1 ! 1=growth limitations, 2=no growth limitations (lim = 1)
	choicePhytoCellKillingScheme          	= 2 ! 1=cell age or sinking below zeu, 2=sinking below 200 m
	
	choiceFractalDimensionScheme          	= 1 ! 1=fixed value, 2=random value
	
	choiceStickinessFunc                  	= 2 ! 1=SLAMS-1.0 (linear function), 2=saturating function with mineral contribution
	choiceMineralEffOnFaecalPellPorosity  	= 1 ! 1=minerals increase porosity of fps, 2=minerals decrease porosity
	choiceCollisionSamplingScheme         	= 2 ! 1=SLAMS-1.0 (no subsampling, stochastic sampling, replacement), 2=SLAMS-2.0 (no subsampling, deterministic sampling, without replacement), 3=SDM (stochastic subsampling, without replacement)
	choiceCollisionKernels                	= 2 ! 1=Burd & Jackson 2009, 2=Jackson 2001
	choiceIncludeHydrodynamicForces       	= 1 ! 1=yes, 0=no (for Burd & Jackson 2009 only!!)
	choiceCriteriaToBreak                 	= 1 ! 1 = based on Kolmogorov length scale, 2 = based on Reynolds number
	
	choiceIsSinkingAffectedByTurbulence   	= .false.
	
	shrinkAfterMineralDissolutionScheme   	= 0 ! 0=no, 1=yes
	shrinkAfterTepDegradationScheme       	= 0 
	shrinkAfterMicrobialMetabolismScheme  	= 0 ! 0=no, 1=yes
 	choiceIsMineralProtectAgainstMicrobResp = .true.
 	
 	choiceZooBehaviourKernels              	= 1 ! 1=unique behaviour, 2=cruisers vs. ambushers
 	choiceZooIngestionCriteria             	= 1 ! 1=volumetric fraction OM, 2=volumetric fraction OM + carbon density + porosity

	! ------------------------------------------------------------------------------------
	! READ NAMELIST FROM FILE
	! ------------------------------------------------------------------------------------

	call findunit(ioUnit) 
	open (unit = ioUnit, file = 'namelist.input', status = 'old', action='read')
		
		! Reference parameters
		rewind(ioUnit)
		read (unit=ioUnit, nml=ReferenceParameters, iostat=ioStatusRef)
		if (ioStatusRef /= 0) then
			write(*,*) 'ERROR reading namelist: ReferenceParameters'
			write(*,*) 'IOSTAT = ', ioStatusRef
			stop
		end if
		
		! Fixed parameters
		rewind(ioUnit)
		read (unit=ioUnit, nml=FixedParameters, iostat=ioStatusFix)
		if (ioStatusFix /= 0) then
			write(*,*) 'ERROR reading namelist: FixedParameters'
			write(*,*) 'IOSTAT = ', ioStatusFix
			stop
		end if
		
		! Configuration parameters
		rewind(ioUnit)
		read (unit=ioUnit, nml=ConfigurationParameters, iostat=ioStatusConf)
		if (ioStatusConf /= 0) then
			write(*,*) 'ERROR reading namelist: ConfigurationParameters'
			write(*,*) 'IOSTAT = ', ioStatusConf
			stop
		end if
		
	close(ioUnit)

	! ------------------------------------------------------------------------------------
	! DERIVED RUNTIME STATE
    !
    ! Compute quantities that depend on runtime configuration and global constants
	! ------------------------------------------------------------------------------------
	
	nTimeStepsDay                      = INT(24d0*3600d0/timeStep) ! dividing the num seconds in a day by nTimeStepsDay has to result in an integer number	
	nTimeStepsYear                     = nTimeStepsDay*365
	maxNumTimeSteps                    = nTimeStepsYear*nYears
	
	nTimeStepsPassedToShiftClusters    = nTimeStepsDay*nDaysPassedToShiftClusters
	nTimeStepsPassedToSeedClay         = nTimeStepsDay*nDaysPassedToSeedClays
	nTimeStepsPassedToDestroyLargeAggs = nTimeStepsDay*nDaysPassedToDestroyLargeAggs
	
	maxNumClustersPerProfile = nNewPhytoClusters + nNewClayClusters + &
							   nNewPhytoTepClusters + maxNumZooDeadClustersPerProfile
	
	C_quota_pft_max = [C_quota_diat_max,   &
					   C_quota_flagel_max, & 
					   C_quota_cocco_max,  &
					   C_quota_pico_max]
						   
	C_quota_pft_min = [C_quota_diat_min,   &
					   C_quota_flagel_min, & 
					   C_quota_cocco_min,  &
					   C_quota_pico_min]
													   
	growth_rate_max_0deg_pft = [growth_rate_max_0deg_diat,   &
					   			growth_rate_max_0deg_flagel, &
					   			growth_rate_max_0deg_cocco,  &
					   			growth_rate_max_0deg_pico] 

	k_NO3_pft = [k_NO3_diat, k_NO3_flagel, k_NO3_cocco, k_NO3_pico]
	
	k_PO4_pft = [k_PO4_diat, k_PO4_flagel, k_PO4_cocco, k_PO4_pico]

	alpha_chl_spec_pft = [alpha_chl_spec_diat,   &
						  alpha_chl_spec_flagel, &
						  alpha_chl_spec_cocco,  &
						  alpha_chl_spec_pico]
						  
	q10_phyto = [q10_diat, q10_flagel, q10_cocco, q10_pico]
		
end subroutine InitialiseRuntimeParameters	

! ========================================================================================

end module modelparameters