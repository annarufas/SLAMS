#include "blockdefinitions.h"

module sink

! ----------------------------------------------------------------------------------------
! This module simulates the gravitational sinking of particles, tracking their movement 
! as they descend through the water column and collecting and imaging them. It also handles 
! the photolysis of particles that accumulate at the surface.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use modelcounters, only: nClustersPhotolysed, nClustersAtSeafloor
use modelconstants, only: SECONDS_PER_DAY, MOLAR_MASS_CARBON, RHO_TEP
use modelparameters, only: timeStep, iPhotoTepC, iCalcite, iOpal, iClay, nPfts, maxNumSmsTerms, &
	maxNumAuxTerms, nVolumeClasses, particleVolumeClasses, choiceIsSinkingAffectedByTurbulence, &
	sea_surface_microlayer_depth, photodegradation_rate_tepc, detection_limit_poc, &
	operational_size_poc_min, shrinkAfterTepDegradationScheme, C_frac_in_TEP, &
	detection_limit_calc, detection_limit_opal, detection_limit_clay, iMicrobSolubOrgC, &
	iMicrobSolubTepC, iMicrobSolubCaCO3, iMicrobSolubOpal, iMicrobSolubClay, choicePhytoCellKillingScheme
use modelgrid, only: gridCellArea
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
use modelparticlecollection, only: lossTerms
#endif
use sanitychecks, only:  WriteStatusAndStop, IsQuantityEffectivelyZero, CheckParticleSanity, &
	IsExceedingInitialAmount
use findfunctions, only: FindParticleAttributeClass
use calcparticleattributes, only: ParticleFractalDimension, ParticleDryMass, ParticleMaterialVolume, &
	ParticleRadius, ParticlePorosity, ParticleDensity, ParticleStickiness, ParticleSettlingVelocity, &
	ShrinkParticle
use primaryproduction, only: KillLivingPhytoplanktonCells
	
implicit none
private
public :: SinkParticle, PhotodegradationOfTep

contains

! ========================================================================================

subroutine SinkParticle(particle, nClusters, nActiveClusters, activeClusterIndices, nSedTrapDeployDepths, &
	sedTrap, sedTrapSeafloor, sedTrapDeployDepths, zeu, localSeafloorDepth, MLD, nLocalDepthLayers, &
	zbot, ztop, muPftProfile, doDiagnosticsAccum, iTimeStep, iProfile)
	
	integer, intent(in) :: nClusters, nActiveClusters, nSedTrapDeployDepths, nLocalDepthLayers, &
		iTimeStep, iProfile
	integer, dimension(nActiveClusters), intent(in) :: activeClusterIndices	
	real*8, intent(in) :: zeu, localSeafloorDepth, MLD
	real*8, dimension(nSedTrapDeployDepths), intent(in) :: sedTrapDeployDepths
	real*8, dimension(nLocalDepthLayers), intent(in) :: zbot, ztop
	real*8, dimension(nPfts,nLocalDepthLayers), intent(in) :: muPftProfile
	logical, intent(in) :: doDiagnosticsAccum
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(5,nSedTrapDeployDepths), intent(inout) :: sedTrap
	real*8, dimension(5,1), intent(inout) :: sedTrapSeafloor

	integer :: iCluster, iActiveCluster, iSedTrapDeployDepth, iVolClass
	real*8 :: particleOldDepth, particleNewDepth, harvest, timeStepDays
	
	if (choiceIsSinkingAffectedByTurbulence) call random_number(harvest)	
	timeStepDays = timeStep/SECONDS_PER_DAY ! time step in days
		
	do iActiveCluster = 1, nActiveClusters
		iCluster = activeClusterIndices(iActiveCluster)
		if (particle(iCluster)%phase /= 1) cycle
		
		! ------ Compute the new particle depth ------
		if (.not. ieee_is_finite(particle(iCluster)%velocity)) then
			write(*,*) 'ERROR: non-finite sinking velocity detected in SinkParticle'
    		call WriteStatusAndStop()
		end if
			
		particleOldDepth = particle(iCluster)%depth
		particleNewDepth = particleOldDepth + (particle(iCluster)%velocity*timeStepDays) ! m
		
		! ------ If particle would be above mixed layer, re-entrain randomly within MLD ------
		if (choiceIsSinkingAffectedByTurbulence) then	
			if (particleNewDepth <= MLD) particleNewDepth = harvest*MLD
		end if

		if (.not. ieee_is_finite(particleNewDepth)) then
			write(*,*) 'ERROR: particle depth is NaN after sinking is applied.'
			write(*,*) '  Old depth:', particleOldDepth
			write(*,*) '  New depth:', particleNewDepth
			write(*,*) '  Sinking velocity:', particle(iCluster)%velocity
			write(*,*) '  Radius:', particle(iCluster)%radius
			write(*,*) '  Density:', particle(iCluster)%density
			write(*,*) '  Org matter content:', particle(iCluster)%massOrgMatter
			write(*,*) '  Porosity:', particle(iCluster)%porosity
			write(*,*) '  Fractal dimension:', particle(iCluster)%fracDim
			call WriteStatusAndStop( )
		end if
					
#ifdef BLOCK_DIAGNOSTICS

	! 	Deploy the sediment trap
	!	Notice that the particle will be recorded at all the depth horizons that it 
	!	has crossed since the last time step. Only positively sinking particles are trapped	
		if (doDiagnosticsAccum) then	
			do iSedTrapDeployDepth = 1, nSedTrapDeployDepths	
				if (particleOldDepth < sedTrapDeployDepths(iSedTrapDeployDepth) .and. &
				  	particleNewDepth > sedTrapDeployDepths(iSedTrapDeployDepth)) then
					call SedimentTrap(nSedTrapDeployDepths, sedTrap(:,:), particle, &
						nClusters, iCluster, iSedTrapDeployDepth, iProfile)
				end if		  	
			end do
		end if
#endif

		! ------- Update particle depth	-------
		particle(iCluster)%depth = particleNewDepth

		! ------ Kill living cells if they sink below zeu or exceed maximum age ------
		if (choicePhytoCellKillingScheme == 1) then ! use zeu and cell age
			call KillLivingPhytoplanktonCells(particle, nClusters, iCluster, zeu, muPftProfile, &
				nLocalDepthLayers, zbot, ztop, iTimeStep, timeStepDays)
		elseif (choicePhytoCellKillingScheme == 2) then ! use the 200 m threshold	
			if (particle(iCluster)%living == 1 .and. particle(iCluster)%depth > 200d0) then
				particle(iCluster)%living = 0
				particle(iCluster)%initType = 4
			end if
		end if
		
		! ------- Bring buoyant particles back to the surface of the water -------
		if (particle(iCluster)%depth < 0d0) particle(iCluster)%depth = 0d0 ! m

		! ------- The particle has reached the seafloor	-------		
		if (particle(iCluster)%depth >= localSeafloorDepth) then
			particle(iCluster)%depth = localSeafloorDepth			 
#ifdef BLOCK_DIAGNOSTICS
            if (doDiagnosticsAccum) then
				  call SedimentTrap(1, sedTrapSeafloor(:,:), particle, nClusters, iCluster, 1, iProfile)
			end if
#endif
			particle(iCluster)%phase = 2							
			nClustersAtSeafloor(iProfile) = nClustersAtSeafloor(iProfile) + 1
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS	
			! Classify particles into volume sections as in Jackson's model				
			iVolClass = FindParticleAttributeClass(nVolumeClasses, particleVolumeClasses, particle(iCluster)%solidVolume)
			if (iVolClass > 0) then
				lossTerms(1,iVolClass) = lossTerms(1,iVolClass)	+ particle(iCluster)%nPxC
			end if
#endif							
		end if ! seafloor
	end do

end subroutine SinkParticle

! ========================================================================================

subroutine SedimentTrap(nDeployDepths, device, particle, nClusters, iCluster, iDeployDepth, iProfile) 

	! ------------------------------------------------------------------------------------
	! Collects the amount of material leaving specific depth horizons (i.e., device 
	! deployment depths). The trap is deployed at every time step from the SinkParticle 
	! subroutine above.
	! ------------------------------------------------------------------------------------

	integer, intent(in) :: nDeployDepths
	real*8, dimension(5,nDeployDepths), intent(inout) :: device	
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	integer, intent(in) :: nClusters, iCluster, iDeployDepth, iProfile	
	
	real*8 :: gramsToMgm2, nParticles

	gramsToMgm2 = 1d3/gridCellArea(iProfile)
	nParticles = particle(iCluster)%nPxC

	device(1,iDeployDepth) = device(1,iDeployDepth) &
		+ particle(iCluster)%molesOrgC*MOLAR_MASS_CARBON*nParticles*gramsToMgm2 ! mg organic carbon m-2 (ATTENTION, what we want here is organic carbon and not organic matter)

	device(2,iDeployDepth) = device(2,iDeployDepth) &
		+ particle(iCluster)%molesTepC*MOLAR_MASS_CARBON*nParticles*gramsToMgm2 ! mg TepC m-2 (ATTENTION, what we want here is TEP carbon and not TEP matter)

	device(3,iDeployDepth) = device(3,iDeployDepth) &
		+ particle(iCluster)%massMineral(iCalcite)*nParticles*gramsToMgm2 ! mg calcite m-2		 

	device(4,iDeployDepth) = device(4,iDeployDepth) &
		+ particle(iCluster)%massMineral(iOpal)*nParticles*gramsToMgm2  ! mg opal m-2

	device(5,iDeployDepth) = device(5,iDeployDepth) &
		+ particle(iCluster)%massMineral(iClay)*nParticles*gramsToMgm2 ! mg clay m-2
		
end subroutine SedimentTrap

! ========================================================================================

subroutine PhotodegradationOfTep(particle, nClusters, SMSterm, nActiveClusters, activeClusterIndices, &
	waterRhoSurface, waterDynViscoSurface, iProfile)

	integer, intent(in) :: nClusters, nActiveClusters, iProfile
	integer, dimension(nActiveClusters), intent(in) :: activeClusterIndices
	real*8, intent(in) :: waterRhoSurface, waterDynViscoSurface
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	
	integer :: iActiveCluster, iCluster
	real*8 :: decayFraction, particleTepC, photoTepC, volDegradedTep, solidVolDegraded
	logical :: carbonZero, tepZero, calciteZero, opalZero, clayZero

	! ------ Decay fraction ------
	decayFraction = 1d0 - EXP(-(photodegradation_rate_tepc/SECONDS_PER_DAY) * timeStep)
	if (.not. ieee_is_finite(decayFraction)) then
    	write(*,*) 'ERROR: decayFraction not finite in PhotodegradationOfTep'
    	call WriteStatusAndStop()
	end if
	decayFraction = MAX(0d0, MIN(1d0, decayFraction))

	do iActiveCluster = 1, nActiveClusters
		iCluster = activeClusterIndices(iActiveCluster)
		if (particle(iCluster)%phase /= 1) cycle
        if (particle(iCluster)%depth > sea_surface_microlayer_depth) cycle
        if (particle(iCluster)%density >= waterRhoSurface) cycle
        particleTepC = particle(iCluster)%molesTepC
		if (particleTepC <= 0d0) cycle
		
		! ------ Degrade TEP ------
		photoTepC = particleTepC * decayFraction ! mol	
		particle(iCluster)%molesTepC = particleTepC - photoTepC					
		SMSterm(iPhotoTepC) = SMSterm(iPhotoTepC) + photoTepC*particle(iCluster)%nPxC
		nClustersPhotolysed(iProfile) = nClustersPhotolysed(iProfile) + 1
		
		! ------ Shrink the particle ------
		if (shrinkAfterTepDegradationScheme > 0 .and. particle(iCluster)%nPpxP >= 2d0) then		
			volDegradedTep = (photoTepC*MOLAR_MASS_CARBON/C_frac_in_TEP)/RHO_TEP ! cm3 
			solidVolDegraded = volDegradedTep*1d12 ! um3 
			if (.not. ieee_is_finite(solidVolDegraded) .or. solidVolDegraded < 0d0) then
    			write(*,*) 'ERROR: invalid degraded solid volume in PhotodegradationOfTep'
    			call WriteStatusAndStop()
			end if
			if (solidVolDegraded > 0d0)  call ShrinkParticle(particle, nClusters, iCluster, solidVolDegraded)
		end if 
		
		! ------ Update attributes -------
		call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
		call ParticleDryMass(particle, nClusters, iCluster)
		call ParticleMaterialVolume(particle, nClusters, iCluster)
		call ParticleRadius(particle, nClusters, iCluster)
		call ParticlePorosity(particle, nClusters, iCluster)
		call ParticleDensity(particle, nClusters, iCluster, waterRhoSurface)
		call ParticleStickiness(particle, nClusters, iCluster)
		call ParticleSettlingVelocity(particle, nClusters, iCluster, waterRhoSurface, waterDynViscoSurface)
	 
	 	! ------- Sanity checks ------
	 	call CheckParticleSanity(particle, nClusters, iCluster, 'TEP degradation')
	 	
	 	! Did we degrade more carbon than existed?
		if (IsExceedingInitialAmount(photoTepC, particleTepC, detection_limit_poc)) then
			write(*,*) 'ERROR: more material is photolysed than the one available.'
			write(*,*) '  Starting and photolysed material:', particleTepC, photoTepC
			write(*,*) '  Radius:', particle(iCluster)%radius
			write(*,*) '  Mineral content:', particle(iCluster)%molesMineral(:)
			write(*,*) '  Depth:', particle(iCluster)%depth
			call WriteStatusAndStop( )
		end if
		
		! Has the particle become too small? --> flush to solubilisation pool and mark as empty
		if (2d0*particle(iCluster)%radius < operational_size_poc_min ) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubCaCO3) = SMSterm(iMicrobSolubCaCO3) + particle(iCluster)%molesMineral(iCalcite)*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubOpal) = SMSterm(iMicrobSolubOpal) + particle(iCluster)%molesMineral(iOpal)*particle(iCluster)%nPxC
			SMSterm(iMicrobSolubClay) = SMSterm(iMicrobSolubClay) + particle(iCluster)%molesMineral(iClay)*particle(iCluster)%nPxC
			particle(iCluster)%phase = 3
			cycle ! nothing else to do here
		end if
		
		! Compute total remaining material to derive a tolerance limit
		carbonZero  = IsQuantityEffectivelyZero(particle(iCluster)%molesOrgC, detection_limit_poc)
		tepZero     = IsQuantityEffectivelyZero(particle(iCluster)%molesTepC, detection_limit_poc)
		calciteZero = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iCalcite), detection_limit_calc)
		opalZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iOpal), detection_limit_opal)
		clayZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iClay), detection_limit_clay)
		
		! Check: did the amount of carbon left go beyond the detection limit? --> allocate leftovers to the solubilised pool
		if (carbonZero) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
			particle(iCluster)%molesOrgC = 0d0
		end if
		if (tepZero) then
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*particle(iCluster)%nPxC
			particle(iCluster)%molesTepC = 0d0	
		end if
	
		! Check: has the particle been emptied after organics disappeared? --> to empty phase
		if (carbonZero .and. tepZero .and. calciteZero .and. opalZero .and. clayZero) then 
			particle(iCluster)%phase = 3
		end if

	end do

end subroutine PhotodegradationOfTep

! ========================================================================================

end module sink