#include "blockdefinitions.h"

module wrappers

! ----------------------------------------------------------------------------------------
! This module is concerned with wrapping the information used by key model processes
! and that is used with the Transport Matrix Method.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use modelcounters, only: particleCollectionPeriod, iTimeStepYear
use modelconstants, only: SECONDS_PER_DAY
use modelparameters, only: nPfts, nDepthLayers, timeStep, nTimeStepsPassedToSeedClay, &
	nTimeStepsPassedToDestroyLargeAggs, maxNumSmsTerms, maxNumAuxTerms, nNewPhytoClusters
use modelforcingdata, only: accumAeolClay
use sanitychecks, only: WriteStatusAndStop
use waterphysicsandlight, only: WaterKinematicViscosity
use primaryproduction, only: PhytoSeedingProbability, PhytoGrowthRate, EuphoticLayerDepth
use injectsurfaceparticles, only: FixedProduction, PhytoplanktonProduction, TepProduction, &
	ClayDeposition
use particledynamics, only: CoagulateParticles, BreakUnstableParticles
use mineraldissolution, only: CalciumCarbonateAndOpalDissolution
use heterotrophicmetabolism, only: MesozooplanktonInteraction, ParticleAttachedMicrobialMetabolism
use modeloutput, only: WriteParticleFlux, WriteParticleAverageAttributes, &
	WriteEulerianVariables, WriteAuxiliaryTerms, WriteAndPrintModelStatistics, WriteTestSmoluchowski
use sink, only: SinkParticle, PhotodegradationOfTep
	 
implicit none
private
public :: ProductionOfSurfaceParticles, ParticleProcessing, GravitationalSettling, &
	WritePeriodicOutput

contains

! ========================================================================================

subroutine ProductionOfSurfaceParticles(particle, nClusters, iLastLocus, zeu, &
	nLocalDepthLayers, ztop, zbot, zmid, gridCellArea, SMSterm, auxTerm, auxCount, &
	muPftProfile, NPP, Chla, PAR0, MLD, AeolClay, NO3, PO4, SiOH4, OmegaCalc, TempC, Rho, &
	DynVisco, iTimeStep, iProfile)        	
	
	integer, intent(in) :: nClusters, nLocalDepthLayers, iTimeStep, iProfile
	real*8, intent(in) :: gridCellArea, NPP, Chla, PAR0, MLD, AeolClay
	real*8, dimension(nLocalDepthLayers), intent(in) :: ztop, zbot, zmid, NO3, PO4, SiOH4, &
		OmegaCalc, TempC, Rho, DynVisco
	integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms,nDepthLayers), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms,nDepthLayers), intent(inout) :: auxTerm, auxCount
	real*8, intent(out) :: zeu
	real*8, dimension(nPfts,nLocalDepthLayers), intent(out) :: muPftProfile

	integer :: nDepthLayersToSeedPhyto
	integer, dimension(nNewPhytoClusters) :: newPhytoClustersIdxs, newPhytoClustersPfts
	real*8 :: molesNewOrgCarbon
	real*8, dimension(nPfts) :: probPftMean
	real*8, dimension(nLocalDepthLayers) :: parz
	real*8, dimension(:,:), allocatable :: probPftProfile
	logical, dimension(nLocalDepthLayers) :: validSeedingDepthIdx

#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	if (iTimeStep == 1) then
		call FixedProduction(particle, nClusters, iLastLocus, SMSterm(:,:), iTimeStep)
	end if

#else

	! ------ Initialise output ------
	muPftProfile(:,:) = 0d0
	
	! ------ Initialise internal variables ------
	probPftMean(:) = 0d0
	molesNewOrgCarbon = 0d0
	zeu = 0d0
	parz(:) = 0d0
	nDepthLayersToSeedPhyto = 0
	validSeedingDepthIdx(:) = .false.

	if (NPP > 0d0 .and. PAR0 > 0d0) then
	
		! Determine euphotic mask, PAR profile and zeu
		if (Chla > 0d0 .and. MLD > 0d0) then
			call EuphoticLayerDepth(zeu, parz(:), nDepthLayersToSeedPhyto, validSeedingDepthIdx(:), &
				auxTerm(:,1), auxCount(:,1), Chla, PAR0, MLD, nLocalDepthLayers, zmid(:)) ! --in primaryproduction.F90
		end if
		
		! Fill growth rates array and calculate seeding probabilities
		if (nDepthLayersToSeedPhyto > 0) then
			call PhytoGrowthRate(muPftProfile(:,:), validSeedingDepthIdx(:), nDepthLayersToSeedPhyto, &
				nLocalDepthLayers, parz(:), NO3(:), PO4(:), SiOH4(:), TempC(:)) ! --in primaryproduction.F90

			allocate(probPftProfile(nDepthLayersToSeedPhyto,nPfts))
			probPftProfile(:,:) = 0d0
			call PhytoSeedingProbability(probPftMean(:), probPftProfile(:,:), muPftProfile(:,:), &
				validSeedingDepthIdx(:), nDepthLayersToSeedPhyto, nLocalDepthLayers, auxTerm(:,:), &
				auxCount(:,:)) ! --in primaryproduction.F90
		
#ifdef BLOCK_PHYTOPLANKTON

			! ------ Seed phytoplankton cells ------
			if (SUM(probPftMean(:)) > 0d0) then
				call PhytoplanktonProduction(newPhytoClustersIdxs(:), newPhytoClustersPfts(:), &
					molesNewOrgCarbon, particle, nClusters, iLastLocus, nLocalDepthLayers, &
					nDepthLayersToSeedPhyto, validSeedingDepthIdx(:), ztop(:), zbot(:), &
					gridCellArea, SMSterm(:,:), auxTerm(:,1), auxCount(:,1), probPftMean(:), &
					probPftProfile(:,:), NPP, OmegaCalc(:), TempC(:), Rho(:), DynVisco(:), iTimeStep) ! --in injectsurfaceparticles.F90
				deallocate(probPftProfile)
			end if	

#ifdef BLOCK_TEP
			
			! ------ Seed TEPs ------
			if (molesNewOrgCarbon > 0d0) then
				call TepProduction(particle, nClusters, iLastLocus, nLocalDepthLayers, &
					ztop(:), zbot(:), SMSterm(:,:), newPhytoClustersIdxs(:), newPhytoClustersPfts(:), &
					Rho(:), DynVisco(:), iTimeStep) ! --in injectsurfaceparticles.F90						
			end if

#endif	
#endif							
		end if ! nDepthLayersToSeedPhyto > 0
	end if ! end checking that there is NPP to seed the model with

#ifdef BLOCK_CLAY
	! ------ Seed clays ------
	! Aeolian deposition events do not happen at every time step, but after intervals of 
	! nDaysPassedToSeedClay (10) days. That does not mean that the dust deposition in the 
	! 10th day is just the corresponding one to that day, but that we have been 
	! accumulating the aeolian dust deposition that there has been in the previous 10 days 
	! and put it altogether in just one cluster in the 10th day.

	if (MOD(iTimeStep,nTimeStepsPassedToSeedClay) == 0) then
		if (accumAeolClay(iProfile) > 0d0) then
			call ClayDeposition(particle, nClusters, iLastLocus, gridCellArea, SMSterm(:,1), &
				accumAeolClay(iProfile), Rho(1), DynVisco(1), iTimeStep) ! --in injectsurfaceparticles.F90			
		end if
		accumAeolClay(iProfile) = 0d0
    else
   	    accumAeolClay(iProfile) = accumAeolClay(iProfile) + AeolClay*timeStep ! g clay m-2
	end if
#endif

#endif

end subroutine ProductionOfSurfaceParticles

! ========================================================================================

subroutine ParticleProcessing(particle, nClusters, iLastLocus, SMSterm, auxTerm, auxCount, &
	nLayerClusters, layerClusterIndices, PAR0, OmegaCalc, O2, TempC, Rho, MesoZoo, DynVisco, &
	tkeLayer, dvmDepthNightUpper, dvmDepthNightLower, dvmDepthDayUpper, dvmDepthDayLower, &
	zmidLayer, gridCellVolume, gridCellArea, iTimeStep, iYear, iProfile)

	integer, intent(in) :: nClusters, nLayerClusters, iTimeStep, iYear, iProfile
	integer, dimension(nLayerClusters), intent(in) :: layerClusterIndices
	real*8, intent(in) :: PAR0, OmegaCalc, O2, TempC, Rho, MesoZoo, DynVisco, tkeLayer, &
		dvmDepthNightUpper, dvmDepthNightLower, dvmDepthDayUpper, dvmDepthDayLower, &
		zmidLayer, gridCellVolume, gridCellArea
	integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount

	real*8 :: TempK, waterRho, waterDynVisco, waterKinVisco, waterShearRate, kolmogorovLengthScale
	
! Hydrodynamic variables are defined at the moment
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	TempK = 20d0 + 273.15d0 ! K			
	waterRho = 1.0275d0 ! g cm-3	
	waterKinVisco = 1d-6 ! m2 s-1
	waterDynVisco = 1d4*waterKinVisco * waterRho ! g cm-1 s-1
	waterShearRate = 1d0 ! s-1	
	kolmogorovLengthScale = (waterKinVisco**3/tkeLayer)**0.25 ! m	
#else
	TempK = TempC + 273.15d0 ! K
	waterDynVisco = DynVisco ! g cm-1 s-1
	waterKinVisco = WaterKinematicViscosity(waterDynVisco,Rho) ! m2 s-1
	waterShearRate = (tkeLayer/waterKinVisco)**0.5 ! s-1 	
	kolmogorovLengthScale = (waterKinVisco**3/tkeLayer)**0.25 ! m	
#endif

	! ------ Sanity checks first ------	
	if (nLayerClusters < 0 .or. nLayerClusters > nClusters) then
    	write(*,*) 'ERROR: invalid nLayerClusters in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(gridCellVolume) .or. gridCellVolume <= 0d0) then
    	write(*,*) 'ERROR: invalid gridCellVolume in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(waterDynVisco) .or. waterDynVisco <= 0d0) then
    	write(*,*) 'ERROR: invalid waterDynVisco in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(Rho) .or. Rho <= 0d0) then
    	write(*,*) 'ERROR: invalid Rho in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(kolmogorovLengthScale) .or. kolmogorovLengthScale <= 0d0) then
    	write(*,*) 'ERROR: invalid kolmogorovLengthScale in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(OmegaCalc)) then
    	write(*,*) 'ERROR: OmegaCalc not finite in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(TempC)) then
		write(*,*) 'ERROR: TempC not finite in ParticleProcessing'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(MesoZoo) .or. MesoZoo <= 0d0) then
    	write(*,*) 'ERROR: invalid MesoZoo in ParticleProcessing'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(O2) .or. O2 <= 0d0) then
    	write(*,*) 'ERROR: invalid O2 in ParticleProcessing'
    	call WriteStatusAndStop()
	end if			
			
#ifdef BLOCK_COAGULATION
	if (nLayerClusters >= 2) then		
		call CoagulateParticles(particle, nClusters, nLayerClusters, layerClusterIndices(:), &
			SMSterm(:), auxTerm(:), auxCount(:), TempK, Rho, waterDynVisco, waterShearRate, &
			kolmogorovLengthScale, gridCellVolume, iProfile) ! --in particledynamics.F90	
	end if			
#endif

#ifdef BLOCK_ZOOPLANKTON
if (iYear > 1 .and. MesoZoo > 0d0) then
	call MesozooplanktonInteraction(particle, nClusters, iLastLocus, SMSterm(:), auxTerm(:), &
		auxCount(:), nLayerClusters, layerClusterIndices(:), MesoZoo, O2, TempC, PAR0, &
		dvmDepthNightLower, dvmDepthDayUpper, Rho, waterDynVisco, waterShearRate, &
		kolmogorovLengthScale, zmidLayer, gridCellVolume, iTimeStep, iProfile) ! --in heterotrophicmetabolism.F90	
end if						
#endif

#ifdef BLOCK_MICROBIAL_METABOLISM
if (iYear > 1) then		
	call ParticleAttachedMicrobialMetabolism(particle, nClusters, SMSterm(:), auxTerm(:), &
		auxCount(:), nLayerClusters, layerClusterIndices(:), O2, TempC, Rho, waterDynVisco, &
		iProfile ) ! --in heterotrophicmetabolism.F90			
end if												
#endif

#ifdef BLOCK_ABIOTIC_DISSOLUTION
if (iYear > 1) then
	call CalciumCarbonateAndOpalDissolution(particle, nClusters, SMSterm(:), nLayerClusters, &
		layerClusterIndices(:), OmegaCalc, TempC, Rho, waterDynVisco, iProfile) ! --in abioticdissolution.F90
end if					
#endif

#ifdef BLOCK_ABIOTIC_FRAGMENTATION
if (iYear > 1) then
	if (MOD(iTimeStep,nTimeStepsPassedToDestroyLargeAggs) == 0) then
		call BreakUnstableParticles(particle, nClusters, nLayerClusters, layerClusterIndices(:), &
			SMSterm(:), Rho, waterDynVisco, kolmogorovLengthScale, iProfile) ! --in particledynamics.F90		
	end if
end if
#endif

end subroutine ParticleProcessing

! ========================================================================================

subroutine GravitationalSettling(particle, nClusters, nActiveClusters, activeClusterIndices, &
	SMSterm, nSedTrapDeployDepths, sedTrap, sedTrapSeafloor, sedTrapDeployDepths, zeu, MLD, &
	localSeafloorDepth, nLocalDepthLayers, ztop, zbot, waterRhoSurface, dynViscoSurface, &
	muPftProfile, doDiagnosticsAccum, iTimeStep, iYear, iProfile)

	integer, intent(in) :: nClusters, nActiveClusters, nSedTrapDeployDepths, nLocalDepthLayers, &
		iTimeStep, iYear, iProfile
	integer, dimension(nActiveClusters), intent(in) :: activeClusterIndices
	real*8, intent(in) :: zeu, MLD, localSeafloorDepth, waterRhoSurface, dynViscoSurface
	real*8, dimension(nLocalDepthLayers) :: ztop, zbot
	real*8, dimension(nSedTrapDeployDepths), intent(in) :: sedTrapDeployDepths
	real*8, dimension(nPfts,nLocalDepthLayers), intent(in) :: muPftProfile
	logical, intent(in) :: doDiagnosticsAccum
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(5,nSedTrapDeployDepths), intent(inout) :: sedTrap
	real*8, dimension(5,1), intent(inout) :: sedTrapSeafloor
	
	! ------ Sanity checks first ------
	if (.not. ieee_is_finite(MLD) .or. MLD < 0d0) then
    	write(*,*) 'ERROR: invalid MLD in GravitationalSettling'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(waterRhoSurface) .or. waterRhoSurface <= 0d0) then
    	write(*,*) 'ERROR: invalid surface density in GravitationalSettling'
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(dynViscoSurface) .or. dynViscoSurface <= 0d0) then
    	write(*,*) 'ERROR: invalid surface viscosity in GravitationalSettling'
    	call WriteStatusAndStop()
	end if


#ifdef BLOCK_SINK

	call SinkParticle(particle, nClusters, nActiveClusters, activeClusterIndices(:), &
		nSedTrapDeployDepths, sedTrap(:,:), sedTrapSeafloor(:,:), sedTrapDeployDepths(:), &
		zeu, localSeafloorDepth, MLD, nLocalDepthLayers, zbot(:), ztop(:), muPftProfile(:,:), &
		doDiagnosticsAccum, iTimeStep, iProfile) ! --in sink.F90			
		
#ifdef BLOCK_PHOTOLYSIS
if (iYear > 2) then
	call PhotodegradationOfTep(particle, nClusters, SMSterm(:), nActiveClusters, activeClusterIndices(:), &
		waterRhoSurface, dynViscoSurface, iProfile) ! --in sink.F90
end if
#endif			
#endif

end subroutine GravitationalSettling

! ========================================================================================

subroutine WritePeriodicOutput(particle, nClusters, iLastLocus, localSeafloorDepth, &
    nAveragingTimeSteps, nProfiles, iTimeStep)

	integer, intent(in) :: nClusters, iLastLocus, nProfiles, iTimeStep, nAveragingTimeSteps
	real*8, intent(in) :: localSeafloorDepth
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
				
	particleCollectionPeriod = dble(nAveragingTimeSteps)*timeStep/SECONDS_PER_DAY ! d (global counter)
	
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	call WriteTestSmoluchowski(nProfiles, iTimeStep)
#else
  
  	! The sediment traps are continuously collecting material, which then is averaged by 
  	! the number of days that the sediment trap has been working. After that, the 
  	! sediment trap is emptied.
	call WriteParticleFlux(localSeafloorDepth,nProfiles) ! flux: mg m-2 / days passed

  	! Average values for size and velo categories of density, velocity, size, stickiness, 
  	! porosity, moles of orgC, CaCO3, opal, clay and TEP			
	call WriteParticleAverageAttributes(nAveragingTimeSteps,localSeafloorDepth,nProfiles)
		
	call WriteEulerianVariables(localSeafloorDepth,nProfiles)
	
#endif

  	call WriteAuxiliaryTerms(localSeafloorDepth,nProfiles)
	call WriteAndPrintModelStatistics(particle, nClusters, iLastLocus, iTimeStep, nProfiles)		
  					
end subroutine WritePeriodicOutput

! ========================================================================================

end module wrappers