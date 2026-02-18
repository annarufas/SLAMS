#include "blockdefinitions.h"

module heterotrophicmetabolism

! ----------------------------------------------------------------------------------------
! This module is concerned with the functions that are related to zooplankton and bacteria.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use safemath, only: SafeDivide, SafeFloorNonNegative
use modelcounters, only: nParticlesEvaluatedForGrazing, nZooFragmentationEvents, &
	nZooDeadClustersPerProfile, nZooIngestionEvents, nFaecalPelletClustersProduced, &
	nZooDeadClustersProduced, nMicrobRespiredClusters
use modelconstants, only: PI, SECONDS_PER_DAY, RHO_ORGMATTER, RHO_TEP, RHO_OPAL, &
	RHO_CALCITE, RHO_CLAY, MOLAR_MASS_CARBON, MOLAR_MASS_OPAL, MOLAR_MASS_CACO3, &
	MOLAR_MASS_CLAY, MOLAR_VOLUME_OXYGEN
use modelparameters, only: nDepthLayers, maxNumSmsTerms, maxNumAuxTerms, timeStep, &
	maxNumZooDeadClustersPerProfile, nNewZooDeadClustersPerDepthLayer, nMesoZooSizeClasses, &
	sumZooWetWeightProducts, zooCarbonWeightQuota, zooWetWeightQuota, zooWetWeightQuotaBinWidth, &
	zooProsomeLength, zooRadius, zooSwimmingSpeed, zoo_distrib_slope, detection_radius_factor_mesozoo, &
	gut_passage_time_mesozoo, dissol_rate_calc, zoo_absorption_eff_carbon, zoo_net_growth_eff, &
	mort_rate_mesozoo, agg_to_zoo_size_ratio, frac_OM_zoo_ingestion_surf, frac_OM_zoo_ingestion_deep, &
	q10_mesozoo, q10_microb, resp_rate_poc_max_0deg_mesozoo, resp_rate_poc_max_0deg_microb, &
	resp_rate_tepc_max_0deg_microb, resp_rate_poc_max_mesozoo, resp_rate_poc_max_microb, &
	resp_rate_tepc_max_microb, k_O2_resp, solub_rate_poc, solub_rate_tepc, C_frac_in_OM, &
	C_frac_in_TEP, operational_size_poc_min, detection_limit_poc, shrinkAfterMicrobialMetabolismScheme, &
	choiceZooBehaviourKernels, choiceZooIngestionCriteria, choiceIsMineralProtectAgainstMicrobResp, &
	iCalcite, iOpal, iClay, iZooNumber, iNightDvmUpperBound, iNightDvmLowerBound, &
	iDayDvmUpperBound, iDayDvmLowerBound, iZooSpecRespRate, iMicrobSpecRespRate, iMicrobSolubOrgC, &
	iMicrobSolubTepC, iMicrobSolubCaCO3, iMicrobSolubOpal, iMicrobSolubClay, iZooIngestOrgC, &
	iZooIngestTepC, iZooEgestOrgC, iZooEgestTepC, iZooRespOrgC, iZooRespTepC, iZooExcretOrgC, &
	iZooExcretTepC, iZooDissolCaCO3, iMicrobRespOrgC, iMicrobRespTepC, iZooDeathOrgC, iZooDeadBiomass, &
	iPrimProdOrgC, iNumParticlesEvalEncounter, iNumParticlesEncountered, iNumParticlesFragmentedZoo, &
	iNumParticlesIngestedZoo, iNumParticlesOmittedZoo, carbon_density_threshold_zoo_ingestion, &
	detection_limit_calc, detection_limit_opal, detection_limit_clay
use sanitychecks, only: IsQuantityEffectivelyZero, CheckParticleSanity, WriteStatusAndStop, &
	IsExceedingInitialAmount					  
use montecarlosampling, only: RandProbCase
use findfunctions, only: FindIndexToNearestLayer, FindDepthLayerMidpointDepth
use calcparticleattributes, only: ParticleFractalDimension, ParticleDryMass, &
	ParticleMaterialVolume, ParticleRadius, ParticlePorosity, ParticleDensity, &
	ParticleStickiness, ParticleSettlingVelocity, ComputeParticleAttributes, ShrinkParticle, &
	VolumetricFractionOfOrganicMatter, FaecalPelletPorosity
use mineraldissolution, only: CaCO3dissolutionLimFactor
use particledynamics, only: BreakParticle

implicit none
private
public :: ComputeBianchiDvmBounds, ZooplanktonDielVerticalMigrationDepth, MesozooplanktonInteraction, &
	ParticleAttachedMicrobialMetabolism

contains

! ========================================================================================

subroutine ComputeBianchiDvmBounds(dvmBoundsBianchi, nSurfLayersBianchi, nMesoLayersBianchi, &
	nLocalDepthLayers, zmid, localSeafloorDepth)
	
	!-------------------------------------------------------------------------------------
	! Finds indexArray that define surface and upper-mesopelagic depth bands.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nLocalDepthLayers
	real*8, intent(in) :: localSeafloorDepth
	real*8, dimension(nLocalDepthLayers), intent(in) :: zmid
	integer, intent(out) :: nSurfLayersBianchi, nMesoLayersBianchi
	integer, dimension(4), intent(out) :: dvmBoundsBianchi
	
	real*8, parameter, dimension(4) :: bianchiTargetDepthBounds = (/0d0, 25d0, 150d0, 500d0/) ! m (surface boundaries and mesopelagic boundaries, from Bianchi et al. 2013)
	
	! Initialisation
	dvmBoundsBianchi(:) = -1
	nSurfLayersBianchi  = 0
	nMesoLayersBianchi  = 0
	
	! Surface band (always exists)
	dvmBoundsBianchi(1) = FindIndexToNearestLayer(bianchiTargetDepthBounds(1),nLocalDepthLayers,zmid)
	dvmBoundsBianchi(2) = FindIndexToNearestLayer(bianchiTargetDepthBounds(2),nLocalDepthLayers,zmid)
	
	nSurfLayersBianchi = dvmBoundsBianchi(2) - dvmBoundsBianchi(1) + 1
	if (nSurfLayersBianchi <= 0) then
		write(*,*) 'SANITY FAIL: invalid surface DVM bounds:', dvmBoundsBianchi(1:2)
		call WriteStatusAndStop()
	end if
	
	! Mesopelagic band (only if deep enough)
	if (localSeafloorDepth >= bianchiTargetDepthBounds(4)) then
		dvmBoundsBianchi(3) = FindIndexToNearestLayer(bianchiTargetDepthBounds(3),nLocalDepthLayers,zmid)
		dvmBoundsBianchi(4) = FindIndexToNearestLayer(bianchiTargetDepthBounds(4),nLocalDepthLayers,zmid)
		
		nMesoLayersBianchi = dvmBoundsBianchi(4) - dvmBoundsBianchi(3) + 1
		if (nMesoLayersBianchi <= 0) then
			write(*,*) 'SANITY FAIL: invalid mesopelagic DVM bounds:', dvmBoundsBianchi(3:4)
			call WriteStatusAndStop()
		end if
	end if
	
	! Safety checks
	if (dvmBoundsBianchi(2) < dvmBoundsBianchi(1)) then
		write(*,*) 'SANITY FAIL: invalid surface bounds ordering', dvmBoundsBianchi
		call WriteStatusAndStop()
	end if

	if (nMesoLayersBianchi > 0) then
		if (dvmBoundsBianchi(4) < dvmBoundsBianchi(3)) then
			write(*,*) 'SANITY FAIL: invalid mesopelagic bounds ordering', dvmBoundsBianchi
			call WriteStatusAndStop()
		end if
	end if
		
end subroutine ComputeBianchiDvmBounds

! ========================================================================================

subroutine ZooplanktonDielVerticalMigrationDepth(dvmDepthNightUpper, dvmDepthNightLower, &
	dvmDepthDayUpper, dvmDepthDayLower, dvmBoundsBianchi, nSurfLayersBianchi, nMesoLayersBianchi, &
	SMSterm, auxTerm, auxCount, nLocalDepthLayers, depthsUpperBoundaries, depthsLowerBoundaries, &
	O2, TempC, Chla, MLD, localSeafloorDepth)

	integer, intent(in) :: nSurfLayersBianchi, nMesoLayersBianchi, nLocalDepthLayers
	integer, dimension(4), intent(in) :: dvmBoundsBianchi
	real*8, intent(in) :: Chla, MLD, localSeafloorDepth
	real*8, dimension(nLocalDepthLayers), intent(in) :: depthsUpperBoundaries, depthsLowerBoundaries, O2, TempC
	real*8, dimension(maxNumSmsTerms,nDepthLayers), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	real*8, intent(out) :: dvmDepthNightUpper, dvmDepthNightLower, dvmDepthDayUpper, dvmDepthDayLower

	integer :: iDepthLayerMaxFood
	real*8 :: avgOxySurf, avgOxyUpperMeso, avgTempSurf, avgTempUpperMeso, gradientOxy, &
		gradientTempC, dayDvmDepth, maxFood, nightDvmDepth, harvestnight, harvestday, &
		nightDvmDepthRange, dayDvmDepthRange
	real*8, dimension(nDepthLayers) :: food
	
	! Day DVM depth, i.e., depth where mesozooplankton migrate downwards at dawn to escape 
	! surface ocean predation
	if (Chla > 0d0 .and. localSeafloorDepth >= 500d0) then
		avgOxySurf       = SUM(O2(dvmBoundsBianchi(1):dvmBoundsBianchi(2)))/dble(nSurfLayersBianchi)
		avgOxyUpperMeso  = SUM(O2(dvmBoundsBianchi(3):dvmBoundsBianchi(4)))/dble(nMesoLayersBianchi)	
		avgTempSurf      = SUM(TempC(dvmBoundsBianchi(1):dvmBoundsBianchi(2)))/dble(nSurfLayersBianchi)
		avgTempUpperMeso = SUM(TempC(dvmBoundsBianchi(3):dvmBoundsBianchi(4)))/dble(nMesoLayersBianchi)	

		gradientOxy = ABS(avgOxySurf-avgOxyUpperMeso)*1d3/MOLAR_VOLUME_OXYGEN ! mL L-1 --> mmol m-3
		gradientTempC = ABS(avgTempSurf-avgTempUpperMeso)

 		! Empirical formula from Bianchi et al. 2013, a function of O2, T, chl a and MLD
		dayDvmDepth = 398d0 - (0.56d0*gradientOxy) - (115d0*LOG10(Chla)) + (0.36d0*MLD) - (2.4d0*gradientTempC) ! m
		if (.not. ieee_is_finite(dayDvmDepth)) then
			write(*,*) 'SANITY FAIL: detected wrong day DVM depth.'
			write(*,*) '  dayDvmDepth:', dayDvmDepth
			write(*,*) '  Chla:', Chla
			write(*,*) '  MLD:', MLD
			write(*,*) '  gradientOxy:', gradientOxy
			write(*,*) '  gradientTempC:', gradientTempC
			call WriteStatusAndStop( )
		end if
		
		! Clamp to [0, seafloor]
        dayDvmDepth = MAX(0d0, MIN(dayDvmDepth, depthsLowerBoundaries(nLocalDepthLayers)))       	
	else
		dayDvmDepth = localSeafloorDepth*0.5d0 ! m
	end if

	! Night DVM depth, i.e., depth where zooplankton ascend when the night falls and coincides
	! with the depth of max food
	food(:) = SMSterm(iPrimProdOrgC,:) 
	maxFood = MAXVAL(food(:))
	if (maxFood <= 0d0) then
        ! no food signal: set night depth near surface (or some default)
        nightDvmDepth = 0d0
    else
		iDepthLayerMaxFood = MAXLOC(food, 1) ! returns index of max in 1-D array
		if (iDepthLayerMaxFood < 1 .or. iDepthLayerMaxFood > nLocalDepthLayers) then
            write(*,*) 'SANITY FAIL: idxMaxFood out of range', iDepthLayerMaxFood
            call WriteStatusAndStop( )
        end if
		nightDvmDepth = MIN(dayDvmDepth, FindDepthLayerMidpointDepth(iDepthLayerMaxFood, &
			nLocalDepthLayers, depthsUpperBoundaries(:), depthsLowerBoundaries(:)))
	end if

	! Random layer thickness between 10 and 100 m
	call random_number(harvestnight)
    nightDvmDepthRange = 10d0 + (100d0-10d0) * harvestnight
    dvmDepthNightUpper = MAX(0d0, nightDvmDepth - 0.5d0*nightDvmDepthRange)
    dvmDepthNightLower = MIN(depthsLowerBoundaries(nLocalDepthLayers), nightDvmDepth + 0.5d0*nightDvmDepthRange)

    call random_number(harvestday)
    dayDvmDepthRange = 10d0 + (100d0-10d0) * harvestday
    dvmDepthDayUpper = MAX(0d0, dayDvmDepth - 0.5d0*dayDvmDepthRange)
    dvmDepthDayLower = MIN(depthsLowerBoundaries(nLocalDepthLayers), dayDvmDepth + 0.5d0*dayDvmDepthRange)

	! Record the depths
	auxTerm(iNightDvmUpperBound) = auxTerm(iNightDvmUpperBound) + dvmDepthNightUpper
	auxTerm(iNightDvmLowerBound) = auxTerm(iNightDvmLowerBound) + dvmDepthNightLower
	auxTerm(iDayDvmUpperBound) 	 = auxTerm(iDayDvmUpperBound) + dvmDepthDayUpper
	auxTerm(iDayDvmLowerBound)   = auxTerm(iDayDvmLowerBound) + dvmDepthDayLower
	
	auxCount(iNightDvmUpperBound) = auxCount(iNightDvmUpperBound) + 1d0
	auxCount(iNightDvmLowerBound) = auxCount(iNightDvmLowerBound) + 1d0
	auxCount(iDayDvmUpperBound)   = auxCount(iDayDvmUpperBound) + 1d0
	auxCount(iDayDvmLowerBound)   = auxCount(iDayDvmLowerBound) + 1d0	
	
end subroutine ZooplanktonDielVerticalMigrationDepth

! ========================================================================================

subroutine MesozooplanktonInteraction(particle, nClusters, iLastLocus, SMSterm, auxTerm, &
	auxCount, nLayerClusters, layerClusterIndices, MesoZoo, O2, TempC, PAR0, dvmDepthNightLower, &
	dvmDepthDayUpper, Rho, waterDynVisco, waterShearRate, kolmogorovLengthScale, zmidLayer, &
	gridCellVolume, iTimeStep, iProfile)
	
	integer, intent(in) :: nClusters, nLayerClusters, iTimeStep, iProfile
	integer, dimension(nLayerClusters), intent(in) :: layerClusterIndices
	real*8, intent(in) :: MesoZoo, O2, TempC, PAR0, dvmDepthNightLower, dvmDepthDayUpper, &
		Rho, waterDynVisco, waterShearRate, kolmogorovLengthScale, zmidLayer, gridCellVolume
	integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount

	integer :: iLayerCluster, iCluster, iZc, iRandZoo, nClustersEncountered, chosenGrazingCategory
	real*8 :: harvest, somaZoo, assimilatedCarbon, startingAmountCarbon, totalZooWetWeight, &
		zooNumbersTotal, zooDistribIntercept, randZooProsomeLength, randZooCarbonQuota, &
		particleRadius, sizeRatio, nParticlesPerCluster, nExpEncOneParticleWithAllZoo, &
		probEncOneParticleAtLeastOnce, probEncAllParticles, fragScore, ingestScore, omitScore, &
		particleEffRadius, orgMatterVolumeFrac, particleCarbonMass, particleBulkVolume, &
		particleCarbonDensity, nExpIngestEvents, nExpFragEvents, nExpOmitEvents, &
		totalExpectedEvents, probChooseIngest, probChooseFrag
	real*8, dimension(nMesoZooSizeClasses) :: zooNumbersByZc, zooDetectionRadiusByZc, &
		reactionDistanceByZc, behaviourKernelByZc, fluidShearKernelByZc, encounterKernelByZc, &
		countZooClassesSelected, probEncPerPair, rateEncPerPair, nExpEncPerPair, &
		ingestScoreByZc, fragScoreByZc, omitScoreByZc, ingestWeightByZc, fragWeightByZc, &
		omitWeightByZc, probIngestByZc
	real*8, parameter :: LOG_UNDERFLOW = log(tiny(1d0)) ! log of the smallest representable REAL*8 number (threshold below which exp(x) underflows to zero in double precision)
	logical :: fragCond, smallEnoughToConsiderIng, ingestCond

	nClustersEncountered = 0 ! bookkeeping variable
	
	zooNumbersByZc(:) 			= 0d0 
	zooDetectionRadiusByZc(:) 	= 0d0
	reactionDistanceByZc(:) 	= 0d0
	behaviourKernelByZc(:) 		= 0d0
	fluidShearKernelByZc(:) 	= 0d0
	encounterKernelByZc(:) 		= 0d0
	probEncPerPair(:)           = 0d0
	rateEncPerPair(:) 		    = 0d0
	nExpEncPerPair(:)           = 0d0
	ingestScoreByZc(:)          = 0d0
	fragScoreByZc(:)            = 0d0
	omitScoreByZc(:)            = 0d0
	ingestWeightByZc(:)         = 0d0
	fragWeightByZc(:)           = 0d0
	omitWeightByZc(:)           = 0d0 
	probIngestByZc(:)           = 0d0
	countZooClassesSelected(:)  = 1d0
	
	! ------ Calculation of the number of zooplankton individuals per m3 ------
	totalZooWetWeight = 10d0**((LOG10(MesoZoo)-(-0.93d0))/0.95d0) ! mg WW m-3 (Kiorboe 2013, Table 2)
	zooDistribIntercept = totalZooWetWeight/sumZooWetWeightProducts ! Z0
	zooNumbersByZc(:) = zooWetWeightQuotaBinWidth(:)*zooDistribIntercept*zooWetWeightQuota(:)**(-zoo_distrib_slope-1d0) ! ind. m-3
	zooNumbersTotal = SUM(zooNumbersByZc(:)) ! ind. m-3	
	if (zooNumbersTotal > 0d0) then
		auxTerm(iZooNumber) = auxTerm(iZooNumber) + zooNumbersTotal ! ind. m-3
		auxCount(iZooNumber) = auxCount(iZooNumber) + 1		
	end if
	
	! ------ Detection radius ------
	zooDetectionRadiusByZc(:) = detection_radius_factor_mesozoo*zooRadius(:) ! um
	
	! ------ Loop through the particles  in the water layer ------
	somaZoo = 0d0 ! mol C

	do iLayerCluster = 1, nLayerClusters
		iCluster = layerClusterIndices(iLayerCluster)
		if (particle(iCluster)%phase /= 1 .or. 2d0*particle(iCluster)%radius < operational_size_poc_min) cycle
		
		particleRadius = particle(iCluster)%radius
		nParticlesPerCluster = particle(iCluster)%nPxC
		if (.not. ieee_is_finite(particleRadius) .or. particleRadius <= 0d0) then
			write(*,*) 'ERROR: invalid radius in MesozooplanktonInteraction'
			call WriteStatusAndStop()
		end if
		if (.not. ieee_is_finite(nParticlesPerCluster) .or. nParticlesPerCluster <= 0d0) then
			write(*,*) 'ERROR: invalid nParticlesPerCluster in MesozooplanktonInteraction', nParticlesPerCluster
			write(*,*) ' nPxP:', particle(iCluster)%nPpxP
			write(*,*) ' radius:', particleRadius
			write(*,*) ' phase:', particle(iCluster)%phase
			call WriteStatusAndStop()
		end if
	
		! Bookkeeping
		startingAmountCarbon = (particle(iCluster)%molesOrgC + particle(iCluster)%molesTepC)*nParticlesPerCluster
		nParticlesEvaluatedForGrazing(iProfile) = nParticlesEvaluatedForGrazing(iProfile) + 1
		auxTerm(iNumParticlesEvalEncounter) = auxTerm(iNumParticlesEvalEncounter) + nParticlesPerCluster
		auxCount(iNumParticlesEvalEncounter) = auxCount(iNumParticlesEvalEncounter) + 1

		! --------- Behaviour kernel ---------
		reactionDistanceByZc(:) = (zooDetectionRadiusByZc(:) + particle(iCluster)%radius)*1d-6 ! m
	
		select case (choiceZooBehaviourKernels)
		case (1) ! unique behaviour kernel
			behaviourKernelByZc(:) = PI * (reactionDistanceByZc(:)**2) * &
				( ((ABS(particle(iCluster)%velocity)/SECONDS_PER_DAY)**2 + 3d0*zooSwimmingSpeed(:)**2) &
				/(3d0*zooSwimmingSpeed(:)) ) ! m3 s-1
	
		case (2) 

			! Fragmenters, or ambushers/passive swimmers (suspended in the water column).
			! Found in surface waters during the day (high prey, high predators) and in 
			! deeper waters during the night
			if ((PAR0 > 0d0 .and. zmidLayer < dvmDepthDayUpper) .or. &
			   (PAR0 == 0d0 .and. zmidLayer >= dvmDepthNightLower)) then
				
				behaviourKernelByZc(:) = PI * (reactionDistanceByZc(:)**2) * &
					(ABS(particle(iCluster)%velocity))/SECONDS_PER_DAY ! m3 s-1
				
			! Consumers/flux feeders, or cruisers/active swimmers. Found in deeper waters 
			! during the day (low prey, low predators) and surface waters during the night
			! (the risk of predation has decreased)							
			else 				
				behaviourKernelByZc(:) = PI * (reactionDistanceByZc(:)**2) * &
					( ((ABS(particle(iCluster)%velocity)/SECONDS_PER_DAY)**2 + 3d0*zooSwimmingSpeed(:)**2) &
					/(3d0*zooSwimmingSpeed(:)) ) ! m3 s-1
			end if
	
		end select
		
		do iZc = 1, nMesoZooSizeClasses			
			if (kolmogorovLengthScale > reactionDistanceByZc(iZc)) then ! laminar flow 
				fluidShearKernelByZc(iZc) = (4d0/3d0) * waterShearRate * (reactionDistanceByZc(iZc)**3) ! m3 s-1 ind-1	
			else ! turbulent flow
				fluidShearKernelByZc(iZc) = 1.3d0 * waterShearRate * (reactionDistanceByZc(iZc)**3) ! m3 s-1 ind-1
			end if
		end do

		!--------- Encounter kernel ---------
		! Volume swept clear by a single real particle on its way to encountering a zooplankton individual of class Zc				
		encounterKernelByZc(:) = behaviourKernelByZc(:) + fluidShearKernelByZc(:) ! m3 s-1
	
		!--------- Mean encounters per single particle by zooplankton class (lambda) ---------
		probEncPerPair(:) = encounterKernelByZc(:) * timeStep / gridCellVolume ! dimensionless
		rateEncPerPair(:) = encounterKernelByZc(:) * zooNumbersByZc(:) ! s-1, # zooplankton of size ri encountered per unit time by the particle in question
		nExpEncPerPair(:) = rateEncPerPair(:) * timeStep ! dimensionless, expected # of encounters a single real particle has with all zooplankton of class Zc during one timestep (=lambda)
		nExpEncOneParticleWithAllZoo = SUM(nExpEncPerPair) ! expected # encounters per single particle with all zooplankton
		
		! --------- Probability that one particle in the cluster is encountered at least once ---------	
		if (nExpEncOneParticleWithAllZoo < 1d-6) then
			! small-lambda expansion to avoid cancellation: 1 - e^{-λ} ≈ λ - λ^2/2
			probEncOneParticleAtLeastOnce = nExpEncOneParticleWithAllZoo - 0.5d0 * nExpEncOneParticleWithAllZoo*nExpEncOneParticleWithAllZoo
		else
			probEncOneParticleAtLeastOnce = 1d0 - exp(-nExpEncOneParticleWithAllZoo)
		end if
		probEncOneParticleAtLeastOnce = MAX(0d0, MIN(1d0, probEncOneParticleAtLeastOnce)) ! clamp [0, 1]
		
		! --------- Probability that every particle in the cluster is encountered at least once ---------
		! Encounters are independent from each other (independent per-particle Poisson processes), 
		! i.e., many zooplankton, sparse interactions. Zooplankton individuals are numerous, 
		! each encounter targets a single particle, and different zooplankton act independently 
		! and rarely encounter the same multiple particles in a way that strongly couples 
		! the events. In such a dilute regime independence can be a good approximation
		
		!probEncAllParticles = probEncOneParticleAtLeastOnce^np, but below we provide a nuemrical stable formula considering probEncOneParticleAtLeastOnce is usually small and np large
		if (probEncOneParticleAtLeastOnce <= 0d0) then
			probEncAllParticles = 0d0
		else if (probEncOneParticleAtLeastOnce >= 1d0) then
			probEncAllParticles = 1d0
		else
			probEncAllParticles = nParticlesPerCluster * log(probEncOneParticleAtLeastOnce) ! = probEncOneParticleAtLeastOnce^np (work in log space to avoid underflow)
			if (probEncAllParticles < LOG_UNDERFLOW) then
				probEncAllParticles = 0d0
			else
				probEncAllParticles = EXP(probEncAllParticles)
			end if
		end if
		
		! Clamp to [0,1]
		if (probEncAllParticles < 0d0) probEncAllParticles = 0d0
		if (probEncAllParticles > 1d0) probEncAllParticles = 1d0

		! --------- DECISION: if probEncAllParticles is high enough --> whole-cluster encountered ---------
		call random_number(harvest)
		if (probEncAllParticles > harvest) then
			nClustersEncountered = nClustersEncountered + 1
			
			do iZc = 1, nMesoZooSizeClasses
				fragScore   = 0d0  ! grows when sizeRatio large and particle is compact (low porosity) and has many primaries
				ingestScore = 0d0  ! grows with volumetric organic fraction and with small sizeRatio (easy to engulf)
				omitScore   = 0d0  ! a fallback or when ingestion/fragmentation unlikely
				
				particleEffRadius = particleRadius * (1d0 - particle(iCluster)%porosity)**(1d0/3d0) ! smaller than the geometric radius for porous aggregates
				sizeRatio = particleEffRadius / zooRadius(iZc)  

				! Conditions
				fragCond = (sizeRatio > agg_to_zoo_size_ratio .and. particle(iCluster)%nPpxP > 1d0) ! fragment if particle is clearly >> zooplankton
				smallEnoughToConsiderIng = (sizeRatio <= agg_to_zoo_size_ratio) ! ingestion allowed if effective size small enough
				
				! Scheme-specific ingestion test
				select case (choiceZooIngestionCriteria)
				case (1) ! ingestion only if organic fraction above threshold
					orgMatterVolumeFrac = VolumetricFractionOfOrganicMatter(particle, nClusters, iCluster) ! palatability metric related to how much of the particle is eatable by zooplankton versus ballast
					if (zmidLayer < 1000d0) then ! waters above 1000 m
						ingestCond = (orgMatterVolumeFrac > frac_OM_zoo_ingestion_surf)
					else ! very deep waters
						ingestCond = (orgMatterVolumeFrac > frac_OM_zoo_ingestion_deep)
					end if
				case (2) ! ingestion only if carbon density above threshold
					particleCarbonMass = (particle(iCluster)%molesOrgC + particle(iCluster)%molesTepC)*MOLAR_MASS_CARBON ! g
					particleBulkVolume = 1d-12*particle(iCluster)%solidVolume / (1d0 - particle(iCluster)%porosity) ! cm3
					particleCarbonDensity = particleCarbonMass / particleBulkVolume ! g cm-3
					ingestCond = (particleCarbonDensity > carbon_density_threshold_zoo_ingestion)
				end select
				
				! Final decision tree for encountered particles
				if (fragCond) then
					fragScore = 1d0
				else if (smallEnoughToConsiderIng) then
					if (ingestCond) then
						ingestScore = 1d0
					else
						omitScore = 1d0
					end if
				else
					omitScore = 1d0
				end if

				ingestScoreByZc(iZc) = ingestScore 
				fragScoreByZc(iZc)   = fragScore   
				omitScoreByZc(iZc)   = omitScore   
				
				! Combine across zooplankton class (per iZc) using class encounter expectation
				ingestWeightByZc(iZc) = nParticlesPerCluster * nExpEncPerPair(iZc) * ingestScoreByZc(iZc)
				fragWeightByZc(iZc)   = nParticlesPerCluster * nExpEncPerPair(iZc) * fragScoreByZc(iZc)
				omitWeightByZc(iZc)   = nParticlesPerCluster * nExpEncPerPair(iZc) * omitScoreByZc(iZc)

			end do
		
			! Expected encounters of outcome-type for entire cluster from class iZc:
			nExpIngestEvents = SUM(ingestWeightByZc(:))
			nExpFragEvents   = SUM(fragWeightByZc(:))
			nExpOmitEvents   = SUM(omitWeightByZc(:))
			totalExpectedEvents = nExpIngestEvents + nExpFragEvents + nExpOmitEvents
			if (.not. ieee_is_finite(totalExpectedEvents) .or. totalExpectedEvents < 0d0) then
				write(*,*) 'ERROR: invalid totalExpectedEvents in MesozooplanktonInteraction'
				call WriteStatusAndStop( )
			end if
		
			! ------ Decide fate by weighted draw ------
			! Derive probabilities (normalise) and sample outcome (fall back to "omit" if no weight) 
			if (totalExpectedEvents <= 0d0) then
				chosenGrazingCategory = 3   ! omit all
			else
				probChooseIngest = nExpIngestEvents / totalExpectedEvents
				probChooseFrag   = nExpFragEvents / totalExpectedEvents
				if (harvest < probChooseIngest) then
					chosenGrazingCategory = 1 ! ingest all
				else if (harvest < (probChooseIngest + probChooseFrag)) then
					chosenGrazingCategory = 2 ! fragment all
				else
					chosenGrazingCategory = 3 ! omit all
				end if
			end if
			
			select case (chosenGrazingCategory)
			case (1) ! ingest all
			
				! Choose a zooplankton size class to attribute ingestion events. For that,
				! derive probabilities (normalise) and choose a size stochastically
				if (nExpIngestEvents > 0d0) then
					probIngestByZc(:) = ingestWeightByZc(:) / nExpIngestEvents
				else
					probIngestByZc(:) = zooNumbersByZc(:) / zooNumbersTotal ! fallback: choose size class proportionally to zooNumberByZc or uniformly
				end if
				iRandZoo = RandProbCase(nMesoZooSizeClasses, probIngestByZc(:))
				randZooProsomeLength = zooProsomeLength(iRandZoo)*1d-3 ! um --> mm
				countZooClassesSelected(iRandZoo) = countZooClassesSelected(iRandZoo) + 1d0
	
				call MesozooplanktonIngestion(assimilatedCarbon, particle, nClusters, &
					SMSterm, auxTerm, auxCount, iCluster, O2, TempC, Rho, waterDynVisco, &
					randZooProsomeLength, iTimeStep, iProfile)

				if (assimilatedCarbon > startingAmountCarbon) then
					write(*,*) 'ERROR: Assimilated C by zooplankton > starting amount.'
					write(*,*) '  Starting and assimilated material:', startingAmountCarbon, assimilatedCarbon
					call WriteStatusAndStop( )
				end if
				
				! This is later used to help select the zooplankton dead body class
				somaZoo = somaZoo + assimilatedCarbon ! mol C
					
				! Log ingested mass							
				SMSterm(iZooIngestOrgC) = SMSterm(iZooIngestOrgC) + particle(iCluster)%molesOrgC*nParticlesPerCluster									
				SMSterm(iZooIngestTepC) = SMSterm(iZooIngestTepC) + particle(iCluster)%molesTepC*nParticlesPerCluster
				auxTerm(iNumParticlesIngestedZoo) = auxTerm(iNumParticlesIngestedZoo) + nParticlesPerCluster
				auxCount(iNumParticlesIngestedZoo) = auxCount(iNumParticlesIngestedZoo) + 1
				
				! Log encounter
				auxTerm(iNumParticlesEncountered) = auxTerm(iNumParticlesEncountered) + nParticlesPerCluster
				auxCount(iNumParticlesEncountered) = auxCount(iNumParticlesEncountered) + 1
			
			case (2) ! frag-all
			
				call BreakParticle(particle, nClusters, iCluster, SMSterm, Rho, waterDynVisco)
				
				! Log fragmented mass
				nZooFragmentationEvents(iProfile) = nZooFragmentationEvents(iProfile) + 1
				auxTerm(iNumParticlesFragmentedZoo) = auxTerm(iNumParticlesFragmentedZoo) + nParticlesPerCluster
				auxCount(iNumParticlesFragmentedZoo) = auxCount(iNumParticlesFragmentedZoo) + 1
				
				! Log encounter
				auxTerm(iNumParticlesEncountered) = auxTerm(iNumParticlesEncountered) + nParticlesPerCluster
				auxCount(iNumParticlesEncountered) = auxCount(iNumParticlesEncountered) + 1

			case (3) ! omit-all
			
				! Log omitted mass
				auxTerm(iNumParticlesOmittedZoo) = auxTerm(iNumParticlesOmittedZoo) + nParticlesPerCluster
				auxCount(iNumParticlesOmittedZoo) = auxCount(iNumParticlesOmittedZoo) + 1
				
				! Log encounter
				auxTerm(iNumParticlesEncountered) = auxTerm(iNumParticlesEncountered) + nParticlesPerCluster
				auxCount(iNumParticlesEncountered) = auxCount(iNumParticlesEncountered) + 1
				
			end select
			
		end if ! probEncAllParticles > harvest
	end do ! loop across clusters in layer
	
	! ------ Form mesozooplankton dead bodies ------
	! Only make dead bodies when there has been ingestion activity and some material ha been assimilated 

	if (somaZoo > 0d0 .and. nZooDeadClustersPerProfile(iProfile) <= maxNumZooDeadClustersPerProfile) then	
		
		! Pick a class probabilistically proportional to counts - if you want the class 
		! most often chosen historically but still allow randomness (useful to avoid 
		! repeatedly selecting the same class when counts are noisy):
		iRandZoo = RandProbCase(nMesoZooSizeClasses, countZooClassesSelected(:))
		randZooCarbonQuota = zooCarbonWeightQuota(iRandZoo)/(1d3*MOLAR_MASS_CARBON) ! mg C ind-1 --> mol C ind-1
	
		call MesozooplanktonDeath(particle, nClusters, SMSterm, auxTerm, auxCount, &
			iLastLocus, MesoZoo, randZooCarbonQuota, Rho, waterDynVisco, zmidLayer, &
			gridCellVolume, iTimeStep, iProfile)

	end if
	
end subroutine MesozooplanktonInteraction

! ========================================================================================

subroutine MesozooplanktonIngestion(clusterAssimilatedCarbon, particle, nClusters, SMSterm, auxTerm, &
	auxCount, iCluster, O2, TempC, Rho, waterDynVisco, randZooProsomeLength, iTimeStep, iProfile)

	! ----------------------------------------------------------------------------------------
	! This will lead to the production of faecal pellets at this specific time step (instantly).
	! If the particle is deemed appetitive and ingested, three main processes take place, in 
	! this order: dissolution of calcium carbonate, packaging of unabsorbed material into 
	! faecal pellets and respiration. The absorbed food goes to satisfy metabolic demands 
	! (growth and reproduction) - gross growth production. There is a fraction that is respired 
	! and what remains is what is used to increase zooplankton tissue - net growth production.
	! ----------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster, iTimeStep, iProfile
	real*8, intent(in) :: O2, TempC, Rho, waterDynVisco, randZooProsomeLength
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	real*8, intent(out) :: clusterAssimilatedCarbon
	
	real*8 :: nParticlesPerCluster, particleOrgC, particleTepC, particleCalc, clusterIngestedCarbon, &
		clusterAbsorbedCarbon, clusterEgestedCarbon, clusterRespiredCarbon, particleRespiredOrgC,  &
		particleRespiredTepC, particleRespiredCarbon, particleAbsorbedOrgC, particleAbsorbedTepC, &
		particleEgestedOrgC, particleEgestedTepC, pressureLimCalc, dissCalc, zooGutCO3ionConc, &
		zooGutVolume, faecalPelletBulkVol, zooOmegaCalc, particlePotentialExcretOrgC, &
		particlePotentialExcretTepC, particlePotentialExcretCarbon, particleExcretOrgC, &
		particleExcretTepC, clusterExcretedCarbon
		
	nParticlesPerCluster = particle(iCluster)%nPxC
	particleOrgC = particle(iCluster)%molesOrgC
	particleTepC = particle(iCluster)%molesTepC  	
	particleCalc = particle(iCluster)%molesMineral(iCalcite) 
	clusterIngestedCarbon = (particleOrgC+particleTepC)*nParticlesPerCluster

	if (.not. ieee_is_finite(particleOrgC) .or. particleOrgC < 0d0 .or. &
    	.not. ieee_is_finite(particleTepC) .or. particleTepC < 0d0) then
    	write(*,*) 'ERROR: invalid particle carbon before ingestion'
    	write(*,*) '  OrgC, TepC:', particleOrgC, particleTepC
    	call WriteStatusAndStop()
	end if

	! ------ Calculate bulk organic matter absorbed and egested ------
	! The absorbed portion of food ingested is that which is absorbed through the gut wall 
	! and taken up into the animal's tissues, with the rest of the food egested as feces. 	
	
	particleAbsorbedOrgC = zoo_absorption_eff_carbon * particleOrgC
	particleAbsorbedTepC = zoo_absorption_eff_carbon * particleTepC
	clusterAbsorbedCarbon = (particleAbsorbedOrgC+particleAbsorbedTepC)*nParticlesPerCluster
	
	particleEgestedOrgC = (1d0 - zoo_absorption_eff_carbon) * particleOrgC
	particleEgestedTepC = (1d0 - zoo_absorption_eff_carbon) * particleTepC
	clusterEgestedCarbon = (particleEgestedOrgC+particleEgestedTepC)*nParticlesPerCluster
	
	particlePotentialExcretOrgC = (1d0 - zoo_net_growth_eff) * particleAbsorbedOrgC
	particlePotentialExcretTepC = (1d0 - zoo_net_growth_eff) * particleAbsorbedTepC
	particlePotentialExcretCarbon = particlePotentialExcretOrgC + particlePotentialExcretTepC	

	! ------ Is respiration possible? ------	
	! Before any further calculations, let's see whether respiration can happen (i.e., 
	! organisms can live in the current environment conditions). If there cannot be 
	! respiration, there won't be metabolic activity of any type.

	particle(iCluster)%molesOrgC = particleAbsorbedOrgC
	particle(iCluster)%molesTepC = particleAbsorbedTepC
					
	call MesozooplanktonRespiration(particleRespiredOrgC, particleRespiredTepC, particle, &
		nClusters, SMSterm, auxTerm, auxCount, iCluster, O2, TempC) ! this function does not update particle attributes		

	particleRespiredCarbon = particleRespiredOrgC + particleRespiredTepC
	clusterRespiredCarbon = particleRespiredCarbon * nParticlesPerCluster

	if (.not. ieee_is_finite(particleRespiredCarbon) .or. particleRespiredCarbon < 0d0) then
    	write(*,*) 'ERROR: invalid particleRespiredCarbon'
    	call WriteStatusAndStop()
	end if
	if (particleRespiredCarbon > (particleAbsorbedOrgC + particleAbsorbedTepC)) then
   		write(*,*) 'ERROR: respiration exceeds absorbed carbon'
    	write(*,*) '  respired, absorbed:', particleRespiredCarbon, particleAbsorbedOrgC+particleAbsorbedTepC
    	call WriteStatusAndStop()
	end if
	
	! ------ Proceed with digestion ------
	if (particleRespiredCarbon > 0d0) then
		nZooIngestionEvents(iProfile) = nZooIngestionEvents(iProfile) + 1
		
		! ------ Calculate volume of the faecal pellet ------
		! From Stamiezskin et al. 2015
		faecalPelletBulkVol = 10d0**(2.58d0*LOG10(randZooProsomeLength) + 5.4d0) ! um3 (PL in mm)
		
		! ------ CaCO3 dissolution in the acidic gut/food vacuole of mesozooplankton ------	
		! Model from Penry & Jumars (1986) and Jumars & Penry (1989), used in Jansen & Wolf-Gladrow (2001)
		if (particleCalc > 0d0) then
		
			! Calculate gut volume of zooplankton (Table 1 of Jansen & Wolf-Gladrow (2001))
			zooGutVolume = 5d0*faecalPelletBulkVol*1d-18 ! um3 --> m3, expecting values around 1-8 x 1e-9 L
			
			! Calculate calcite dissolution limitation factor
			if (zooGutVolume <= 0d0) then
			  	zooOmegaCalc = 0d0
			else
			  	zooGutCO3ionConc = particleCalc*1d3 / zooGutVolume ! mmol m-3
			  	zooOmegaCalc = zooGutCO3ionConc / 42.7d0 ! 42.7 mmol m-3 (=umol kg-1) is a surface ocean value (100 m), as in Jansen's thesis
			end if
			pressureLimCalc = CaCO3dissolutionLimFactor(zooOmegaCalc)
			
			! Dissolve calcite		
			dissCalc = particleCalc &
				- particleCalc*EXP(- ( (dissol_rate_calc*pressureLimCalc/SECONDS_PER_DAY) &
					+ (1d0/gut_passage_time_mesozoo) ) * timeStep) ! mol
			particle(iCluster)%molesMineral(iCalcite) = particleCalc - dissCalc
			SMSterm(iZooDissolCaCO3) = SMSterm(iZooDissolCaCO3) + dissCalc*nParticlesPerCluster
		end if

		! ------ Excretory losses (urine) ------	
		if (clusterIngestedCarbon > 0d0 .and. particleRespiredCarbon < particlePotentialExcretCarbon) then
			particleExcretOrgC = particlePotentialExcretOrgC - particleRespiredOrgC
			particleExcretTepC = particlePotentialExcretTepC - particleRespiredTepC
			clusterExcretedCarbon = (particleExcretOrgC+particleExcretTepC) * nParticlesPerCluster
			SMSterm(iZooExcretOrgC) = SMSterm(iZooExcretOrgC) + particleExcretOrgC*nParticlesPerCluster
			SMSterm(iZooExcretTepC) = SMSterm(iZooExcretTepC) + particleExcretTepC*nParticlesPerCluster			
		else		
			clusterExcretedCarbon = 0d0				
		end if

		! ------ Faecal pellet production ------
		! Use the minerals that are left and the material that has not been absorbed
		particle(iCluster)%molesOrgC = particleEgestedOrgC
		particle(iCluster)%molesTepC = particleEgestedTepC
		if ((particleEgestedOrgC+particleEgestedTepC+SUM(particle(iCluster)%molesMineral(:))) > 0d0) then	
			call MesozooplanktonEgestion(particle, nClusters, iCluster, SMSterm, &
				faecalPelletBulkVol, Rho, waterDynVisco, iTimeStep, iProfile) ! this function updates particle attributes
		end if
		
		SMSterm(iZooEgestOrgC) = SMSterm(iZooEgestOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
		SMSterm(iZooEgestTepC) = SMSterm(iZooEgestTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
							
		! ------ Assimilate leftover material ------
		clusterAssimilatedCarbon = clusterAbsorbedCarbon - clusterRespiredCarbon - clusterExcretedCarbon  ! mol C
		if (.not. ieee_is_finite(clusterAssimilatedCarbon) .or. clusterAssimilatedCarbon < 0d0) then
    		write(*,*) 'ERROR: invalid clusterAssimilatedCarbon', clusterAssimilatedCarbon
    		call WriteStatusAndStop()
		end if
		if (clusterAssimilatedCarbon > clusterAbsorbedCarbon) then
			write(*,*) 'ERROR: assimilated carbon exceeds absorbed'
			call WriteStatusAndStop()
		end if

	else

		particle(iCluster)%molesOrgC = particleOrgC
		particle(iCluster)%molesTepC = particleTepC
		clusterAssimilatedCarbon = 0d0

	end if ! end checking whether metabolism can be active in current temperature and O2 conditions

end subroutine MesozooplanktonIngestion

! ========================================================================================

subroutine MesozooplanktonRespiration(particleRespiredOrgC, particleRespiredTepC, particle, &
	nClusters, SMSterm, auxTerm, auxCount, iCluster, O2, TempC)

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: O2, TempC
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	real*8, intent(out) :: particleRespiredOrgC, particleRespiredTepC
	
	real*8 :: nParticlesPerCluster, particleOrgC, particleTepC, tempFunc, oxyLim, &
		respRateCarbon, startingCarbon, particleRespiredCarbon

	! ------ Initialise output ------
	particleRespiredOrgC = 0d0
	particleRespiredTepC = 0d0
	
	! ------ Sanity checks ------
	nParticlesPerCluster = particle(iCluster)%nPxC
	particleOrgC = particle(iCluster)%molesOrgC
	particleTepC = particle(iCluster)%molesTepC
	startingCarbon = particleOrgC + particleTepC
	if (.not. ieee_is_finite(particleOrgC) .or. particleOrgC < 0d0 .or. &
    	.not. ieee_is_finite(particleTepC) .or. particleTepC < 0d0) then
    	write(*,*) 'ERROR: invalid particle carbon before MesozooplanktonRespiration'
    	write(*,*) '  OrgC, TepC:', particleOrgC, particleTepC
    	call WriteStatusAndStop()
	end if
	
	! ------ Temperature and oxygen scaling ------
    tempFunc = q10_mesozoo**((TempC-0d0)/10d0)
	oxyLim = O2/(O2+k_O2_resp)
	respRateCarbon = MIN( (resp_rate_poc_max_0deg_mesozoo*tempFunc/SECONDS_PER_DAY), &
		(resp_rate_poc_max_mesozoo/SECONDS_PER_DAY) ) * oxyLim ! s-1
	if (.not. ieee_is_finite(respRateCarbon) .or. respRateCarbon < 0d0) then
    	write(*,*) 'ERROR: invalid respRateCarbon in MesozooplanktonRespiration', respRateCarbon
    	call WriteStatusAndStop()
	end if
	
	! ------ Save diagnostics ------
	auxTerm(iZooSpecRespRate) = auxTerm(iZooSpecRespRate) + respRateCarbon
	auxCount(iZooSpecRespRate) = auxCount(iZooSpecRespRate) + 1	
	
	! ------ Respire ------
	if (particleOrgC > 0d0) then		
		particleRespiredOrgC = particleOrgC - particleOrgC*EXP(-respRateCarbon*timeStep) ! mol C (in the exponential solution, respOrgC cannot be > particleOrgC)
		particle(iCluster)%molesOrgC = particleOrgC - particleRespiredOrgC
		SMSterm(iZooRespOrgC) = SMSterm(iZooRespOrgC) + particleRespiredOrgC*nParticlesPerCluster
	end if
	if (particleTepC > 0d0) then
		particleRespiredTepC = particleTepC - particleTepC*EXP(-respRateCarbon*timeStep) ! mol C (in the exponential solution, respTepC cannot be > particleTepC)
		particle(iCluster)%molesTepC = particleTepC - particleRespiredTepC
		SMSterm(iZooRespTepC) = SMSterm(iZooRespTepC) + particleRespiredTepC*nParticlesPerCluster
	end if
	particleRespiredCarbon = particleRespiredOrgC + particleRespiredTepC
	
	! ------ Sanity checks ------
	if (particle(iCluster)%molesOrgC < 0d0 .or. particle(iCluster)%molesTepC < 0d0) then
    	write(*,*) 'ERROR: negative particle carbon after MesozooplanktonRespiration'
    	call WriteStatusAndStop()
	end if
	if (IsExceedingInitialAmount(particleRespiredCarbon, startingCarbon, detection_limit_poc)) then
		write(*,*) 'ERROR: more material is respired by mesozooplankton than the one available.'
		write(*,*) '  Starting and respired material:', startingCarbon, particleRespiredCarbon
		write(*,*) '  tempFunc and oxyLim', tempFunc, oxyLim
		write(*,*) '  Radius:', particle(iCluster)%radius
		write(*,*) '  Mineral content:', particle(iCluster)%molesMineral(:)
		write(*,*) '  Depth:', particle(iCluster)%depth	
    	call WriteStatusAndStop( )
	end if

end subroutine MesozooplanktonRespiration

! ========================================================================================

subroutine MesozooplanktonEgestion(particle, nClusters, iCluster, SMSterm, faecalPelletBulkVol, &
	Rho, waterDynVisco, iTimeStep, iProfile)

	integer, intent(in) :: nClusters, iCluster, iTimeStep, iProfile
	real*8, intent(in) :: faecalPelletBulkVol, Rho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm	
	
	integer :: safetyCounter
	real*8 :: nPpxParentParticle, nPxParentCluster, parentMolesOrgC, parentMolesTepC, &
		parentMolesOpal, parentMolesCalc, parentMolesClay, parentClusterMass, parentSolidVol, &
		parentSolidVolOrgMatter, parentSolidVolTep, parentSolidVolOpal, parentSolidVolCalc, &
		parentSolidVolClay, faecalPorosity, faecalSolidVol, faecalToParentSolidVolume, &
		faecalSolidVolOrgMatter, faecalSolidVolTep, faecalSolidVolOpal, faecalSolidVolCalc, &
		faecalSolidVolClay, faecalClusterMass, pelletMass, expectedNumPellets, realisedNumPellets, &
		expectedNumPrimParticles, realisedNumPrimParticles, nPrimParticles, remainder, harvest, &
		calcFaecalPelletSolidVolFromMass
	logical :: isFaecalIdenticalToParent
	
	! ------ Bookkeeping ------		
	nFaecalPelletClustersProduced(iProfile) = nFaecalPelletClustersProduced(iProfile) + 1

	! ------ Update these ------
	particle(iCluster)%faecal = 1 ! the particle is turned into faecal material
	particle(iCluster)%living = 0 ! the particle is doomed dead
	particle(iCluster)%initType = 1 ! faecal type
	particle(iCluster)%tstepCreat = iTimeStep
	particle(iCluster)%phase = 1			
	call ParticleDryMass(particle, nClusters, iCluster)
	call ParticleMaterialVolume(particle, nClusters, iCluster)
	
	! Check: is the food particle about to be processed for egestion physically consistent?
	call CheckParticleSanity(particle, nClusters, iCluster, 'just before faecal pellet creation')

	! ------ Record the attributes of the parent particle ------	
	nPpxParentParticle = particle(iCluster)%nPpxP
	nPxParentCluster = particle(iCluster)%nPxC
	parentMolesOrgC = particle(iCluster)%molesOrgC
	parentMolesTepC = particle(iCluster)%molesTepC
	parentMolesOpal = particle(iCluster)%molesMineral(iOpal)
	parentMolesCalc = particle(iCluster)%molesMineral(iCalcite)
	parentMolesClay = particle(iCluster)%molesMineral(iClay)
	parentClusterMass = particle(iCluster)%mass*nPxParentCluster			
	parentSolidVol = particle(iCluster)%solidVolume	! um3					
	parentSolidVolOrgMatter = particle(iCluster)%massOrgMatter/RHO_ORGMATTER ! cm3
	parentSolidVolTep = particle(iCluster)%massTep/RHO_TEP ! cm3
	parentSolidVolOpal = particle(iCluster)%massMineral(iOpal)/RHO_OPAL ! cm3
	parentSolidVolCalc = particle(iCluster)%massMineral(iCalcite)/RHO_CALCITE ! cm3
	parentSolidVolClay = particle(iCluster)%massMineral(iClay)/RHO_CLAY ! cm3
	if (.not. ieee_is_finite(parentSolidVol) .or. parentSolidVol <= 0d0) then
    	write(*,*) 'ERROR: invalid parentSolidVol in MesozooplanktonEgestion', parentSolidVol
    	call WriteStatusAndStop()
	end if
	
	! ------ Scaling factor ------
	! The solid volume of the faecal particle is used to scale the material content of the 
	! ingested particle to that of the faecal particle (we use the solid volume and not, for 
	! instance, the radius because the radius implies already a porosity, which might be very 
	! different in the ingested and faecal particles). The faecal particle and the parent 
	! particle are compared in terms of volume. Use the ratio to compute moles, solid volume 
	! and mass of the faecal particles.
	faecalPorosity = FaecalPelletPorosity(particle, nClusters, iCluster)
	faecalSolidVol = (1d0-faecalPorosity)*faecalPelletBulkVol ! um3
	faecalToParentSolidVolume = faecalSolidVol/parentSolidVol
	if (.not. ieee_is_finite(faecalSolidVol) .or. faecalSolidVol <= 0d0) then
    	write(*,*) 'ERROR: invalid faecalSolidVol in MesozooplanktonEgestion', faecalSolidVol
    	call WriteStatusAndStop()
	end if

	! Scale down (or up) the content of the specific materials that are in the parent particles
	faecalSolidVolOrgMatter = parentSolidVolOrgMatter * faecalToParentSolidVolume ! cm3
	faecalSolidVolTep = parentSolidVolTep * faecalToParentSolidVolume ! cm3
	faecalSolidVolOpal = parentSolidVolOpal * faecalToParentSolidVolume ! cm3
	faecalSolidVolCalc = parentSolidVolCalc * faecalToParentSolidVolume ! cm3
	faecalSolidVolClay = parentSolidVolClay * faecalToParentSolidVolume ! cm3

	particle(iCluster)%molesOrgC = faecalSolidVolOrgMatter*RHO_ORGMATTER*C_frac_in_OM/MOLAR_MASS_CARBON
	particle(iCluster)%molesTepC = faecalSolidVolTep*RHO_TEP*C_frac_in_TEP/MOLAR_MASS_CARBON
	particle(iCluster)%molesMineral(iOpal) = faecalSolidVolOpal*RHO_OPAL/MOLAR_MASS_OPAL
	particle(iCluster)%molesMineral(iCalcite) = faecalSolidVolCalc*RHO_CALCITE/MOLAR_MASS_CACO3
	particle(iCluster)%molesMineral(iClay) = faecalSolidVolClay*RHO_CLAY/MOLAR_MASS_CLAY

	! ------ Mass of the faecal pellet ------
	! Calculate the mass of the particle and, from the relationship mass cluster parent to 
	! mass faecal particle, we can work out how many faecal particles per cluster we will create.	
	call ParticleDryMass(particle, nClusters, iCluster)
	pelletMass = particle(iCluster)%mass
	if (.not. ieee_is_finite(pelletMass) .or. pelletMass <= 0d0) then
		write(*,*) 'ERROR: invalid pelletMass', pelletMass
		call WriteStatusAndStop( )
	end if
	
	isFaecalIdenticalToParent = .false.
	
! ----------------------------------------------------------------------------------------
! Scenario A: pelletMass > parentClusterMass --> pellet is too large for parent cluster 
! --> repeatedly shrink the pellet until mass matches (i.e., faecalClusterMass = parentClusterMass) 
! or pellet volume <= parent. If pellet becomes smaller than parent, make faecal identical 
! to parent and exit.
! ----------------------------------------------------------------------------------------

	if (pelletMass > parentClusterMass) then 
		particle(iCluster)%nPxC = 1d0
		faecalClusterMass = pelletMass * particle(iCluster)%nPxC
	
		safetyCounter = 0
		do while (ABS(parentClusterMass-faecalClusterMass) > 1d-12 .and. safetyCounter < 200)
			safetyCounter = safetyCounter + 1
			
			! ------ Shrink per-pellet solid volumes by factor 0.1 ------
			! Stop scaling down either when the faecal cluster mass becomes the same as the 
			! parent cluster mass, or when the faecal pellet volume becomes smaller than that 
			! of the parent (in this latter case, we make the faecal particles identical as the
			! parent particles and then we exit the loop)
	
			faecalSolidVolOrgMatter = faecalSolidVolOrgMatter * 1d-1 ! cm3
			faecalSolidVolTep       = faecalSolidVolTep       * 1d-1 ! cm3
			faecalSolidVolOpal      = faecalSolidVolOpal      * 1d-1 ! cm3
			faecalSolidVolCalc      = faecalSolidVolCalc      * 1d-1 ! cm3
			faecalSolidVolClay      = faecalSolidVolClay      * 1d-1 ! cm3
		
			faecalSolidVol = 1d12*(faecalSolidVolOrgMatter + faecalSolidVolTep &
				+ faecalSolidVolOpal + faecalSolidVolCalc + faecalSolidVolClay) ! um3
		
			particle(iCluster)%molesOrgC = faecalSolidVolOrgMatter*RHO_ORGMATTER*C_frac_in_OM/MOLAR_MASS_CARBON
			particle(iCluster)%molesTepC = faecalSolidVolTep*RHO_TEP*C_frac_in_TEP/MOLAR_MASS_CARBON
			particle(iCluster)%molesMineral(iOpal) = faecalSolidVolOpal*RHO_OPAL/MOLAR_MASS_OPAL
			particle(iCluster)%molesMineral(iCalcite) = faecalSolidVolCalc*RHO_CALCITE/MOLAR_MASS_CACO3
			particle(iCluster)%molesMineral(iClay) = faecalSolidVolClay*RHO_CLAY/MOLAR_MASS_CLAY
	
			! ------ Recompute num pellets ------
			call ParticleDryMass(particle, nClusters, iCluster)
			pelletMass = particle(iCluster)%mass
			
			expectedNumPellets = SafeDivide(parentClusterMass, pelletMass, 0d0)
			if (.not. ieee_is_finite(expectedNumPellets) .or. expectedNumPellets < 0d0) then
   	 			write(*,*) 'ERROR: invalid expectedNumPellets'
    			write(*,*) ' expectedNumPellets=', expectedNumPellets
    			write(*,*) ' parentClusterMass=', parentClusterMass
    			write(*,*) ' pelletMass=', pelletMass
    			call WriteStatusAndStop()
			end if
			realisedNumPellets = SafeFloorNonNegative(expectedNumPellets) ! round down
			particle(iCluster)%nPxC = MAX(1d0, realisedNumPellets)
			
			if (pelletMass < detection_limit_poc*MOLAR_MASS_CARBON) then
    			write(*,*) 'WARNING: pelletMass below detection limit:', pelletMass
    			write(*,*) '	Pellet making process in Scenario A stops and faecal pellet adopts characteristics of parent material'
    			write(*,*) ' 	faecalSolidVol=', faecalSolidVol
    			write(*,*) ' 	parentSolidVol=', parentSolidVol
    			write(*,*) ' 	expectedNumPellets=', expectedNumPellets
			end if

			! NOTICE: we don't want to round up (which could happen if we used stochastic rounding)
			! as that will incur in extra material we don't have to form primary particles)

			! ------ Re-check ------
			! If, in this cycle of making the faecal particles smaller, the faecal particle 
			! becomes smaller in volume than the parent particle, we will end up making 
			! faecal particles made out of just one (primary) particle --> not ideal. In 
			! this case, make the faecal particle identical to the parent particle and stop 
			! calculations.
							
			if (faecalSolidVol <= parentSolidVol) then			
				particle(iCluster)%molesOrgC = parentMolesOrgC
				particle(iCluster)%molesTepC = parentMolesTepC
				particle(iCluster)%molesMineral(iOpal) = parentMolesOpal
				particle(iCluster)%molesMineral(iCalcite) = parentMolesCalc
				particle(iCluster)%molesMineral(iClay) = parentMolesClay 
				particle(iCluster)%nPxC = nPxParentCluster
				particle(iCluster)%nPpxP = nPpxParentParticle
				particle(iCluster)%nPpxC = particle(iCluster)%nPpxP * particle(iCluster)%nPxC

				call ParticleDryMass(particle, nClusters, iCluster)
				call ParticleMaterialVolume(particle, nClusters, iCluster)
				call ParticleFractalDimension( particle, nClusters, iCluster, particle(iCluster)%initType)
				call ParticleRadius(particle, nClusters, iCluster)
				call ParticlePorosity(particle, nClusters, iCluster)
				call ParticleDensity(particle, nClusters, iCluster, Rho)
				call ParticleStickiness(particle, nClusters, iCluster)	
				call ParticleSettlingVelocity(particle, nClusters, iCluster, Rho, waterDynVisco)
							
				isFaecalIdenticalToParent = .true.	
				exit
			end if
			
			! ------ Safety: break if pelletMass becomes zero or non-finite ------
			if (.not. ieee_is_finite(pelletMass) .or. pelletMass <= 0d0) then
				write(*,*) 'ERROR: pelletMass invalid while shrinking'
				call WriteStatusAndStop( )
			end if
					
		end do ! end shrink loop
		
		! If we left due to safetyCounter limit, handle gracefully
    	if (safetyCounter >= 200) then
        	write(*,*) 'WARNING: shrink loop reached iteration limit; proceeding with best-effort pellet sizing'
    	end if
		
		! If we have managed to have mass cluster parent = mass cluster faecal, without 
		! the faecal particle volume going too small (vol faecal > vol parent), finalise 
		! particle fields from the last faecalSolidVol
		
		if (.not. isFaecalIdenticalToParent) then

			particle(iCluster)%solidVolume = faecalSolidVol ! um3 
			!call ParticleDryMass( particle, nClusters, iCluster ) ! no need for it, already calculated	
			
			call ParticleMaterialVolume(particle, nClusters, iCluster)
			calcFaecalPelletSolidVolFromMass = particle(iCluster)%solidVolume
			particle(iCluster)%solidVolume = faecalSolidVol
			if (faecalSolidVol /= calcFaecalPelletSolidVolFromMass) then
				write(*,*) 'WARNING: faecal pellet solid volume calculated normally and from function differ'
				write(*,*) '  Calculated normally, calculated from function:', calcFaecalPelletSolidVolFromMass, faecalSolidVol
			end if	
			
			call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
			call ParticleStickiness(particle, nClusters, iCluster)
			particle(iCluster)%porosity = FaecalPelletPorosity(particle, nClusters, iCluster) ! different from standard calculation of particle porosity
			if (particle(iCluster)%porosity /= faecalPorosity) then
				write(*,*) 'WARNING: faecal pellet porosity pre-calculated differs from new', faecalPorosity, particle(iCluster)%porosity 
				call WriteStatusAndStop( )
			end if
			
			call ParticleDensity(particle, nClusters, iCluster, Rho) ! g cm-3
			
			! NOTICE: particle radius is calculated after porosity, and not the other 
			! way around (as we do for other particle types) since we have not worked out 
			! yet the number of primary particles per faecal aggregate
			if (.not. ieee_is_finite(particle(iCluster)%solidVolume) .or. particle(iCluster)%solidVolume <= 0d0) then
   	 			write(*,*) 'ERROR: invalid faecal solid volume'
    			call WriteStatusAndStop()
			end if
			if (.not. ieee_is_finite(particle(iCluster)%porosity) .or. &
    			particle(iCluster)%porosity < 0d0 .or. particle(iCluster)%porosity >= 1d0) then
    			write(*,*) 'ERROR: invalid faecal porosity'
    			call WriteStatusAndStop()
			end if
			particle(iCluster)%radius = ((3d0*particle(iCluster)%solidVolume)/(4d0*PI*(1d0-particle(iCluster)%porosity)))**(1d0/3d0) ! um	
			particle(iCluster)%radiusPp = particle(iCluster)%radius / ((1d0-particle(iCluster)%porosity)**(1d0/(particle(iCluster)%fracDim-3d0)))
			
			call ParticleSettlingVelocity(particle, nClusters, iCluster, Rho, waterDynVisco) ! m d-1 
			
			! Number of primary particles in the faecal particle (we can use stochastic rounding safely for primaries)
			expectedNumPrimParticles = (particle(iCluster)%radius/particle(iCluster)%radiusPp)**particle(iCluster)%fracDim
			if (.not. ieee_is_finite(expectedNumPrimParticles) .or. expectedNumPrimParticles < 0d0) then
    			write(*,*) 'ERROR: invalid expectedNumPrimParticles in faecal pellet production', expectedNumPrimParticles
    			call WriteStatusAndStop()
			end if
			
			realisedNumPrimParticles = SafeFloorNonNegative(expectedNumPrimParticles) ! round down
			remainder = expectedNumPrimParticles - realisedNumPrimParticles ! used as the probability of adding one extra particle
			remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
			call random_number(harvest)
			if (harvest < remainder) then
				nPrimParticles = realisedNumPrimParticles + 1d0
			else
				nPrimParticles = realisedNumPrimParticles
			end if
			if (nPrimParticles < 1d0) then
				write(*,*) 'INTERNAL ERROR: nPrimParticles < 1 in faecal pellet production'
				call WriteStatusAndStop()
			end if
	
			particle(iCluster)%nPpxP = nPrimParticles
			particle(iCluster)%nPpxC = particle(iCluster)%nPpxP * particle(iCluster)%nPxC
			
			! Compute faecal cluster mass and the remainder, and send any remaining (physical) 
			! leftover mass to microbial solubilisation
			faecalClusterMass = pelletMass * particle(iCluster)%nPxC
			remainder = parentClusterMass - faecalClusterMass ! g
			if (IsQuantityEffectivelyZero(remainder, detection_limit_poc*MOLAR_MASS_CARBON)) remainder = 0d0
			if (remainder > 0d0) then
				SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + remainder/MOLAR_MASS_CARBON
			end if

		end if
		
! ----------------------------------------------------------------------------------------
! Scenario B: the mass of the cluster parent is bigger than the mass of a faecal particle 
! --> produce integer number of pellets.
! ----------------------------------------------------------------------------------------

	elseif (pelletMass <= parentClusterMass) then
	
		! NOTICE: the following sequence is slightly different to the above one as we need
		! to compute particle(iCluster)%nPxC.
	
		particle(iCluster)%solidVolume = faecalSolidVol ! um3
		particle(iCluster)%porosity = faecalPorosity
		!call ParticleDryMass( particle, nClusters, iCluster ) ! no need for it, already calculated
		call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
		call ParticleStickiness(particle, nClusters, iCluster)
		call ParticleDensity(particle, nClusters, iCluster, Rho) ! g cm-3
		particle(iCluster)%radius = ((3d0*particle(iCluster)%solidVolume)/(4d0*PI*(1d0-particle(iCluster)%porosity)))**(1d0/3d0) ! um
		particle(iCluster)%radiusPp = particle(iCluster)%radius/((1d0-particle(iCluster)%porosity)**(1d0/(particle(iCluster)%fracDim-3d0)))
		call ParticleSettlingVelocity(particle, nClusters, iCluster, Rho, waterDynVisco) ! m d-1

		! Compute integer number of pellets producible by parent cluster mass
		expectedNumPellets = SafeDivide(parentClusterMass, pelletMass, 0d0)
		realisedNumPellets = SafeFloorNonNegative(expectedNumPellets) ! round down
		particle(iCluster)%nPxC = MAX(1d0, realisedNumPellets)
		
		! NOTICE: we don't want to round up (which could happen if we used stochastic rounding)
		! as that will incur in extra material we don't have to form particles)

		! Calculate other particle counts (we can safely used stochastic rounding here)
		expectedNumPrimParticles = (particle(iCluster)%radius/particle(iCluster)%radiusPp)**particle(iCluster)%fracDim
		if (.not. ieee_is_finite(expectedNumPrimParticles) .or. expectedNumPrimParticles < 0d0) then
			write(*,*) 'ERROR: invalid expectedNumPrimParticles in faecal pellet production', expectedNumPrimParticles
			call WriteStatusAndStop()
		end if
		
		realisedNumPrimParticles = SafeFloorNonNegative(expectedNumPrimParticles) ! round down
		remainder = expectedNumPrimParticles - realisedNumPrimParticles ! used as the probability of adding one extra particle
		remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
		call random_number(harvest)
		if (harvest < remainder) then
			nPrimParticles = realisedNumPrimParticles + 1d0
		else
			nPrimParticles = realisedNumPrimParticles
		end if
		if (nPrimParticles < 1d0) then
    		write(*,*) 'ERROR: nPrimParticles < 1 in faecal pellet production'
   	 		call WriteStatusAndStop()
		end if

		particle(iCluster)%nPpxP = nPrimParticles
		particle(iCluster)%nPpxC = particle(iCluster)%nPpxP * particle(iCluster)%nPxC

		! Compute faecal cluster mass and the remainder, and send any remaining (physical) 
		! leftover mass to microbial solubilisation
		faecalClusterMass = pelletMass * particle(iCluster)%nPxC
		remainder = parentClusterMass - faecalClusterMass ! g
		if (IsQuantityEffectivelyZero(remainder, detection_limit_poc*MOLAR_MASS_CARBON)) remainder = 0d0
		if (remainder > 0d0) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + remainder/MOLAR_MASS_CARBON
		end if

	end if
	
	! Check: is the particle physically consistent?
	call CheckParticleSanity(particle, nClusters, iCluster, 'faecal pellet creation')

end subroutine MesozooplanktonEgestion

! ========================================================================================

subroutine MesozooplanktonDeath(particle, nClusters, SMSterm, auxTerm, auxCount, iLastLocus, &
	MesoZoo, randZooCarbonWeightQuota, Rho, waterDynVisco, zmidLayer, gridCellVolume, &
	iTimeStep, iProfile)

	! ------------------------------------------------------------------------------------
    ! Some zooplankton will die due to starvation, hatchling failure, disease, etc. It is a 
    ! prognostic process and a function of the current standing zoo biomass. This subroutine 
    ! creates zooplankton carcasses from organic matter considered "dead". One cluster only
	! is created, which will contain all the carcasses.
    ! ------------------------------------------------------------------------------------
    
    integer, intent(in) :: nClusters,  iTimeStep, iProfile	
    real*8, intent(in) :: MesoZoo, randZooCarbonWeightQuota, Rho, waterDynVisco, zmidLayer, gridCellVolume
    integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	
	integer :: iCluster, iZooCluster
	real*8 :: expectedNumZooDeadBodiesPerCluster, realisedNumZooDeadBodiesPerCluster, &
		nZooDeadBodiesPerCluster, harvest, zooLivingBiomass, temporaryZooDeadBiomass, &
		updatedZooDeadBiomass

	! ------ Compute number of dead carcasses that could be formed ------
	zooLivingBiomass = MesoZoo/(1d3*MOLAR_MASS_CARBON) ! mg C m-3 --> mol C m-3
	
	temporaryZooDeadBiomass = zooLivingBiomass*gridCellVolume &
		- (zooLivingBiomass*gridCellVolume)*EXP((-mort_rate_mesozoo/SECONDS_PER_DAY)*timeStep) ! mol C

	if (.not. ieee_is_finite(temporaryZooDeadBiomass) .or. temporaryZooDeadBiomass < 0d0) then
		write(*,*) 'ERROR: invalid temporaryZooDeadBiomass'
		write(*,*) '  zooLivingBiomass=', zooLivingBiomass
		write(*,*) '  gridCellVolume=', gridCellVolume
		write(*,*) '  mort_rate_mesozoo=', mort_rate_mesozoo
		write(*,*) '  timeStep=', timeStep
		call WriteStatusAndStop()
	end if
	
	if (temporaryZooDeadBiomass < randZooCarbonWeightQuota) return ! not enough material to form a dead body

	expectedNumZooDeadBodiesPerCluster = SafeDivide(temporaryZooDeadBiomass, randZooCarbonWeightQuota*dble(nNewZooDeadClustersPerDepthLayer), 0d0)
	if (.not. ieee_is_finite(expectedNumZooDeadBodiesPerCluster) .or. expectedNumZooDeadBodiesPerCluster < 0d0) then
		write(*,*) 'ERROR: invalid expectedNumZooDeadBodiesPerCluster', expectedNumZooDeadBodiesPerCluster
		write(*,*) '  temporaryZooDeadBiomass=', temporaryZooDeadBiomass
		write(*,*) '  randZooCarbonWeightQuota=', randZooCarbonWeightQuota
		call WriteStatusAndStop()
	end if
	
	realisedNumZooDeadBodiesPerCluster = SafeFloorNonNegative(expectedNumZooDeadBodiesPerCluster) ! round down
	
	! NOTICE: we don't want to round up (which could happen if we used stochastic rounding)
	! as that will incur in extra material we don't have to form particles)
	
	nZooDeadBodiesPerCluster = realisedNumZooDeadBodiesPerCluster
	if (nZooDeadBodiesPerCluster < 1d0) return
	
	! ------ Form dead carcasses ------	
	nZooDeadClustersPerProfile(iProfile) = nZooDeadClustersPerProfile(iProfile) + nNewZooDeadClustersPerDepthLayer	

	do iZooCluster = 1, nNewZooDeadClustersPerDepthLayer	
		nZooDeadClustersProduced(iProfile) = nZooDeadClustersProduced(iProfile) + 1	
			
		! ------ Assign next free slot ------
		if (iLastLocus + 1 > nClusters) then
            write(*,*) 'ERROR: particle array exhausted in MesozooplanktonDeath'
            call WriteStatusAndStop()
        end if
		iCluster = iLastLocus + 1
		
		! ------ Compute the attributes ------ 	
		call ComputeParticleAttributes(particle, nClusters, iCluster, 7, 0, 0, 0, iTimeStep, &
			0d0, zmidLayer, nZooDeadBodiesPerCluster, randZooCarbonWeightQuota, Rho, waterDynVisco)
			
		! ------ Store diagnostics ------	
		SMSterm(iZooDeathOrgC) = SMSterm(iZooDeathOrgC) + particle(iCluster)%molesOrgC*nZooDeadBodiesPerCluster ! mol orgC
		
		! ------ Update on success seeding ------
		iLastLocus = iLastLocus + 1
	end do	
	
	! ------ Update 'dead' pool ------
	updatedZooDeadBiomass = nZooDeadBodiesPerCluster*randZooCarbonWeightQuota*dble(nNewZooDeadClustersPerDepthLayer)
	if (updatedZooDeadBiomass > 0d0) then
		auxTerm(iZooDeadBiomass) = auxTerm(iZooDeadBiomass) + updatedZooDeadBiomass/gridCellVolume ! mol C m-3				
		auxCount(iZooDeadBiomass) = auxCount(iZooDeadBiomass) + 1
	end if
								
	! ------ Sanity check ------
	if (IsExceedingInitialAmount(updatedZooDeadBiomass, temporaryZooDeadBiomass, detection_limit_poc)) then
		write(*,*) 'ERROR: Updated zoo dead biomass > temporary (initially estimated).'
		write(*,*) '  Living biomass', zooLivingBiomass
		write(*,*) '  Zoo C quota', randZooCarbonWeightQuota
		write(*,*) '  Updated dead biomass', updatedZooDeadBiomass
		write(*,*) '  Initial dead biomass',temporaryZooDeadBiomass
		write(*,*) '  nZooDeadBodiesPerCluster', nZooDeadBodiesPerCluster
		call WriteStatusAndStop( )
	end if
	
end subroutine MesozooplanktonDeath

! ========================================================================================

subroutine ParticleAttachedMicrobialMetabolism(particle, nClusters, SMSterm, auxTerm, &
	auxCount, nLayerClusters, layerClusterIndices, O2, TempC, Rho, waterDynVisco, iProfile)

	! ------------------------------------------------------------------------------------
    ! Particle-attached microbial metabolism encompass bacteria and microzooplankton. 
    ! When they deal with a particle, two processes happen: respiration and solubilisation.
    ! ------------------------------------------------------------------------------------
 	
 	integer, intent(in) :: nClusters, nLayerClusters, iProfile
 	integer, dimension(nLayerClusters), intent(in) :: layerClusterIndices
 	real*8, intent(in) :: O2, TempC, Rho, waterDynVisco 
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount

	integer :: iLayerCluster, iCluster
	real*8 :: tempFunc, oxyLim, respRatePoc, respRateTepC, solubFracPoc, solubFracTepC
	
	! ------ Precompute respiration capped rates (s-1) ------
	tempFunc = q10_microb**(TempC/10d0)
	oxyLim = O2/(O2+k_O2_resp)
	
	respRatePoc = ( MIN( resp_rate_poc_max_0deg_microb * tempFunc, resp_rate_poc_max_microb ) / &
                SECONDS_PER_DAY ) * oxyLim
    respRateTepC = ( MIN( resp_rate_tepc_max_0deg_microb * tempFunc, resp_rate_tepc_max_microb ) / &
                SECONDS_PER_DAY ) * oxyLim
                
    if (.not. ieee_is_finite(respRatePoc) .or. respRatePoc < 0d0) then
    	write(*,*) 'ERROR: invalid respRatePoc in ParticleAttachedMicrobialMetabolism'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(respRateTepC) .or. respRateTepC < 0d0) then
    	write(*,*) 'ERROR: invalid respRateTepC in ParticleAttachedMicrobialMetabolism'
    	call WriteStatusAndStop()
	end if

	! ------ Solubilisation fraction over this timestep ------
    solubFracPoc = 1d0 - EXP( -( solub_rate_poc / SECONDS_PER_DAY ) * timeStep )
    solubFracTepC = 1d0 - EXP( -( solub_rate_tepc / SECONDS_PER_DAY ) * timeStep )
    solubFracPoc = MAX(0d0, MIN(1d0, solubFracPoc))
	solubFracTepC = MAX(0d0, MIN(1d0, solubFracTepC))

	do iLayerCluster = 1, nLayerClusters
		iCluster = layerClusterIndices(iLayerCluster)
	
		! ------ Respiration (may change phase or reduce POC) ------
		if ( particle(iCluster)%phase == 1 .and. &
			(particle(iCluster)%molesOrgC > 0d0 .or. particle(iCluster)%molesTepC > 0d0) ) then
			
			call MicrobialRespiration( particle, nClusters, SMSterm, auxTerm, auxCount, iCluster, &
				respRatePoc, respRateTepC, Rho, waterDynVisco, iProfile )
		
		end if
		
		! ------- Solubilisation -------
		! Check status again as it might have changed after respiration
		if ( particle(iCluster)%phase == 1 .and. &
			(particle(iCluster)%molesOrgC > 0d0 .or. particle(iCluster)%molesTepC > 0d0) ) then
			
			call MicrobialSolubilisation( particle, nClusters, SMSterm, auxTerm, auxCount, iCluster, &
				solubFracPoc, solubFracTepC, Rho, waterDynVisco, iProfile )
			
		end if
			
	end do

end subroutine ParticleAttachedMicrobialMetabolism

! ========================================================================================

subroutine MicrobialRespiration(particle, nClusters, SMSterm, auxTerm, auxCount, iCluster, &
	respRatePoc, respRateTepC, Rho, waterDynVisco, iProfile)

	integer, intent(in) :: nClusters, iCluster, iProfile
	real*8, intent(in) :: respRatePoc, respRateTepC, Rho, waterDynVisco 
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	
	real*8 :: labilityFactor, particleOrgC, particleTepC, particleRespiredCarbon, respOrgC, &
		respTepC, volDegradedOrgMatter, volDegradedTep, solidVolDegraded, startingAmountCarbon, &
		orgMatterVolumeFrac
	logical :: carbonZero, tepZero, calciteZero, opalZero, clayZero

	! ------ Diagnostics ------
	! Store only the respiration rate for orgC (leave TEP-C aside)
	auxTerm(iMicrobSpecRespRate) = auxTerm(iMicrobSpecRespRate) + respRatePoc
	auxCount(iMicrobSpecRespRate) = auxCount(iMicrobSpecRespRate) + 1

	! ------ Compute a lability factor (the opposite of a protection factor) ------
	if (choiceIsMineralProtectAgainstMicrobResp) then			
		orgMatterVolumeFrac = VolumetricFractionOfOrganicMatter(particle, nClusters, iCluster)
		labilityFactor = orgMatterVolumeFrac
	else
		labilityFactor = 1d0
	end if
		
	! ------ Respire ------		
	particleOrgC         = particle(iCluster)%molesOrgC
	particleTepC         = particle(iCluster)%molesTepC
	startingAmountCarbon = particleOrgC + particleTepC
	respOrgC             = 0d0
	respTepC             = 0d0
			
	if (particleOrgC > 0d0) then
		respOrgC = particleOrgC - particleOrgC*EXP(-labilityFactor*respRatePoc*timeStep) ! mol (in the exponential solution, respOrgC cannot be > particleOrgC)		
		particleOrgC = particleOrgC - respOrgC					
		SMSterm(iMicrobRespOrgC) = SMSterm(iMicrobRespOrgC) + respOrgC*particle(iCluster)%nPxC	
	end if
	if (particleTepC > 0d0) then
		respTepC = particleTepC - particleTepC*EXP(-labilityFactor*respRateTepC*timeStep) ! mol (in the exponential solution, respTepC cannot be > particleTepC)	
		particleTepC = particleTepC - respTepC						   	
		SMSterm(iMicrobRespTepC) = SMSterm(iMicrobRespTepC) + respTepC*particle(iCluster)%nPxC 
	end if 
	particleRespiredCarbon = respOrgC + respTepC

	if (particleRespiredCarbon > 0d0) then
		nMicrobRespiredClusters(iProfile) = nMicrobRespiredClusters(iProfile) + 1
	
		! Update the organic material
		particle(iCluster)%molesOrgC = particleOrgC				
		particle(iCluster)%molesTepC = particleTepC
		
		! ------ Shrink ------	
		! The shrink function will distribute the remaining material into a reduced
		! number of primary particles (if no shrinking is applied, the primary particle 
		! count does not decrease, the primary particles only get smaller)	
		if (shrinkAfterMicrobialMetabolismScheme > 0 .and. particle(iCluster)%nPpxP >= 2d0) then
			volDegradedOrgMatter = (respOrgC*MOLAR_MASS_CARBON/C_frac_in_OM)/RHO_ORGMATTER ! cm3			
			volDegradedTep = (respTepC*MOLAR_MASS_CARBON/C_frac_in_TEP)/RHO_TEP ! cm3 
			select case (shrinkAfterMicrobialMetabolismScheme)
			case (1) ! shrink if both orgC and TEP-C are degraded
				solidVolDegraded = (volDegradedOrgMatter+volDegradedTep)*1d12 ! um3 
			case (2) ! shrink if only TEP-C is degraded
				solidVolDegraded = volDegradedTep*1d12 ! um3 
			end select
			call ShrinkParticle(particle, nClusters, iCluster, solidVolDegraded)
		end if 	
		
		! ------ Update attributes ------
		call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
		call ParticleDryMass(particle, nClusters, iCluster) 
		call ParticleMaterialVolume(particle, nClusters, iCluster) 
		call ParticleRadius(particle, nClusters, iCluster) 
		call ParticlePorosity(particle, nClusters, iCluster)
		call ParticleDensity(particle, nClusters, iCluster, Rho)
		call ParticleStickiness(particle, nClusters, iCluster)
		call ParticleSettlingVelocity(particle, nClusters, iCluster, Rho, waterDynVisco)
 
		! ------ Sanity checks ------
		
		! Did more carbon than the one available get respired?
		if (IsExceedingInitialAmount(particleRespiredCarbon, startingAmountCarbon, detection_limit_poc)) then
			write(*,*) 'ERROR: more material is respired by microbes than the one available.'
			write(*,*) '  Starting and solub. material:', startingAmountCarbon, particleRespiredCarbon
			write(*,*) '  Radius:', particle(iCluster)%radius
			write(*,*) '  Mineral content:', particle(iCluster)%molesMineral(:)
			write(*,*) '  Depth:', particle(iCluster)%depth
			call WriteStatusAndStop( )
		end if
		
		! Is the particle physically consistent?
		call CheckParticleSanity(particle, nClusters, iCluster, 'microbial respiration') 
		
		! Has the particle become too small? --> flush to solubilisation pool and mark as empty
		if (2d0*particle(iCluster)%radius < operational_size_poc_min ) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubCaCO3) = SMSterm(iMicrobSolubCaCO3) + particle(iCluster)%molesMineral(iCalcite)*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubOpal) = SMSterm(iMicrobSolubOpal) + particle(iCluster)%molesMineral(iOpal)*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubClay) = SMSterm(iMicrobSolubClay) + particle(iCluster)%molesMineral(iClay)*particle(iCluster)%nPxC
			particle(iCluster)%phase = 3
			return ! nothing else to do here
		end if
		
		! Compute total remaining material to derive a tolerance limit
		carbonZero  = IsQuantityEffectivelyZero(particle(iCluster)%molesOrgC, detection_limit_poc)
		tepZero     = IsQuantityEffectivelyZero(particle(iCluster)%molesTepC, detection_limit_poc)
		calciteZero = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iCalcite), detection_limit_calc)
		opalZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iOpal), detection_limit_opal)
		clayZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iClay), detection_limit_clay)
		
		! Did the amount of carbon left go beyond the detection limit? --> allocate leftovers to the solubilised pool
		if (carbonZero) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
			particle(iCluster)%molesOrgC = 0d0
		end if
		if (tepZero) then
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
			particle(iCluster)%molesTepC = 0d0	
		end if

		! Has the particle been emptied after organics disappeared? --> to empty phase
		if (carbonZero .and. tepZero .and. calciteZero .and. opalZero .and. clayZero) then 
			particle(iCluster)%phase = 3
		end if
		
	end if ! particleRespiredCarbon > 0 

end subroutine MicrobialRespiration

! ========================================================================================

subroutine MicrobialSolubilisation(particle, nClusters, SMSterm, auxTerm, auxCount, iCluster, &
	solubFracPoc, solubFracTepC, Rho, waterDynVisco, iProfile)

	integer, intent(in) :: nClusters, iCluster, iProfile
	real*8, intent(in) :: solubFracPoc, solubFracTepC, Rho, waterDynVisco 
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount

	real*8 :: particleOrgC, particleTepC, particleSolubilisedCarbon, solubOrgC, solubTepC, &
		volDegradedOrgMatter, volDegradedTep, solidVolDegraded, startingAmountCarbon
	logical :: carbonZero, tepZero, calciteZero, opalZero, clayZero

	! ------ Solubilise ------
	particleOrgC         = particle(iCluster)%molesOrgC
	particleTepC         = particle(iCluster)%molesTepC
	startingAmountCarbon = particleOrgC + particleTepC
	solubOrgC            = 0d0
	solubTepC            = 0d0
						
	if (particleOrgC > 0d0) then
		solubOrgC = particleOrgC*solubFracPoc
		particleOrgC = particleOrgC - solubOrgC
		SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + solubOrgC*particle(iCluster)%nPxC
	end if
	if (particleTepC > 0d0) then
		solubTepC = particleTepC*solubFracTepC ! mol
		particleTepC = particleTepC - solubTepC
		SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + solubTepC*particle(iCluster)%nPxC
	end if
	particleSolubilisedCarbon = solubOrgC + solubTepC
	
	if (particleSolubilisedCarbon > 0d0) then
	
		! Update the organic material
		particle(iCluster)%molesOrgC = particleOrgC				
		particle(iCluster)%molesTepC = particleTepC	
		
		! ------ Shrink ------
		! The shrink function will distribute the remaining material into a reduced
		! number of primary particles (if no shrinking is applied, the primary particle 
		! count does not decrease, the primary particles only get smaller)	
		if (shrinkAfterMicrobialMetabolismScheme > 0 .and. particle(iCluster)%nPpxP >= 2d0) then
			volDegradedOrgMatter = (solubOrgC*MOLAR_MASS_CARBON/C_frac_in_OM)/RHO_ORGMATTER ! cm3			
			volDegradedTep = (solubTepC*MOLAR_MASS_CARBON/C_frac_in_TEP)/RHO_TEP ! cm3 
			select case (shrinkAfterMicrobialMetabolismScheme)
			case (1) ! shrink if both orgC and TEP-C are degraded
				solidVolDegraded = (volDegradedOrgMatter+volDegradedTep)*1d12 ! um3 
			case (2) ! shrink if only TEP-C is degraded
				solidVolDegraded = volDegradedTep*1d12 ! um3 
			end select
			call ShrinkParticle( particle, nClusters, iCluster, solidVolDegraded )
		end if 	
		
		! ------ Update attributes ------
		call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
		call ParticleDryMass(particle, nClusters, iCluster)
		call ParticleMaterialVolume(particle, nClusters, iCluster) 
		call ParticleRadius(particle, nClusters, iCluster)
		call ParticlePorosity(particle, nClusters, iCluster)
		call ParticleDensity(particle, nClusters, iCluster, Rho)
		call ParticleStickiness(particle, nClusters, iCluster)
		call ParticleSettlingVelocity(particle, nClusters, iCluster, Rho, waterDynVisco)
 
		! ------ Sanity checks ------
		
		! Did more carbon than the one available get solubilised?
		if (IsExceedingInitialAmount(particleSolubilisedCarbon, startingAmountCarbon, detection_limit_poc)) then
			write(*,*) 'ERROR: more material is solubilised by microbes than the one available.'
			write(*,*) '  Starting and solub. material:', startingAmountCarbon, particleSolubilisedCarbon
			write(*,*) '  Radius:', particle(iCluster)%radius
			write(*,*) '  Mineral content:', particle(iCluster)%molesMineral(:)
			write(*,*) '  Depth:', particle(iCluster)%depth
			call WriteStatusAndStop( )
		end if
		
		! Is the particle physically consistent?
		call CheckParticleSanity(particle, nClusters, iCluster, 'microbial solubilisation') 
		
		! Check: has the particle become too small? --> flush to solubilisation pool and mark as empty
		if (2d0*particle(iCluster)%radius < operational_size_poc_min ) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubCaCO3) = SMSterm(iMicrobSolubCaCO3) + particle(iCluster)%molesMineral(iCalcite)*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubOpal) = SMSterm(iMicrobSolubOpal) + particle(iCluster)%molesMineral(iOpal)*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubClay) = SMSterm(iMicrobSolubClay) + particle(iCluster)%molesMineral(iClay)*particle(iCluster)%nPxC
			particle(iCluster)%phase = 3
			return ! nothing else to do here
		end if
		
		! Compute total remaining material to derive a tolerance limit
		carbonZero  = IsQuantityEffectivelyZero(particle(iCluster)%molesOrgC, detection_limit_poc)
		tepZero     = IsQuantityEffectivelyZero(particle(iCluster)%molesTepC, detection_limit_poc)
		calciteZero = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iCalcite), detection_limit_calc)
		opalZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iOpal), detection_limit_opal)
		clayZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iClay), detection_limit_clay)
				
		! Did the amount of carbon left go beyond the detection limit? --> allocate leftovers to the solubilised pool
		if (carbonZero) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
			particle(iCluster)%molesOrgC = 0d0
		end if
		if (tepZero) then
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
			particle(iCluster)%molesTepC = 0d0	
		end if
	
		! Has the particle been emptied after organics disappeared? --> to empty phase
		if (carbonZero .and. tepZero .and. calciteZero .and. opalZero .and. clayZero) then 
			particle(iCluster)%phase = 3
		end if
		
	end if ! particleSolubilisedCarbon > 0

end subroutine MicrobialSolubilisation

! ========================================================================================

end module heterotrophicmetabolism