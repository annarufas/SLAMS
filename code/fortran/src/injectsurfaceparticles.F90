module injectsurfaceparticles

! ----------------------------------------------------------------------------------------
! This module is concerned with the creation of primary particles (i.e., non aggregated
! particles) at the surface ocean: phytoplankton cells, transparent exopolymer particles
! (TEP) and clays.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use safemath, only: SafeNINTNonNegative, SafeDivide, SafeFloorNonNegative
use modelconstants, only: MOLAR_MASS_CARBON, MOLAR_MASS_CLAY
use modelparameters, only: nDepthLayers, nNewPhytoClusters, nNewPhytoTepClusters, &
	nNewClayClusters, nNewFixedClusters, nPfts, maxNumSmsTerms, maxNumAuxTerms, timeStep, &
	choicePhytoSeedingDepthDistrib, C_quota_pft_max, C_quota_pft_min, clay_quota, &
	phyto_exudation_frac, C_quota_TEP_max, C_quota_TEP_min, iPrimProdOrgC, iProdTepPhyto, &
	iPrimProdCaCO3, iPrimProdOpal, iDepoClay, iCalcite, iOpal, iClay, iFreshDiatCellQuota, &
	iFreshFlagelCellQuota, iFreshCoccoCellQuota, iFreshPicoCellQuota, iDiatBiomass, &
	iFlagelBiomass, iCoccoBiomass, iPicoBiomass, iMicrobSolubOrgC, iMicrobSolubClay, &
	iMicrobSolubTepC
use sanitychecks, only: CheckParticleSanity, PrintParticleProperties, WriteStatusAndStop 
use montecarlosampling, only: RandProbCase
use findfunctions, only: FindDepthLayerIndex, FindDepthLayerMidpointDepth
use calcparticleattributes, only: ComputeParticleAttributes, ParticleDryMass, &
	ParticleMaterialVolume, ParticleRadius, ParticleDensity, ParticleSettlingVelocity

implicit none
private
public :: PhytoplanktonProduction, TepProduction, ClayDeposition, FixedProduction

contains

! ========================================================================================

subroutine PhytoplanktonProduction(newPhytoClustersIdxs, newPhytoClustersPfts, molesNewOrgCarbon, &
	particle, nClusters, iLastLocus, nLocalDepthLayers, nDepthLayersToSeedPhyto, validSeedingDepthIdx, &
	ztop, zbot, gridCellArea, SMSterm, auxTerm, auxCount, probPftMean, probPftProfile, NPP, &
	OmegaCalc, TempC, Rho, DynVisco, iTimeStep)
	
	! ------------------------------------------------------------------------------------
	! This subroutine is called every time step that there's light available only. It
	! converts continuous primary production (NPP) into discrete phytoplankton particles while
	! conserving organic carbon in expectation.
	!
	! Net primary production (NPP) is a continuous flux (mol C m-2 s-1), whereas phytoplankton
	! particles are discrete entities. Converting between the two requires enforcing
	! integer particle counts without introducing systematic carbon loss or gain. We enforce 
	! discreteness using stochastic rounding:
	!   - The expected (real-valued) number of particles per cluster is computed.
	!   - The integer number of particles is obtained by rounding up or down with a
	!     probability equal to the fractional remainder.
	!   - Carbon is therefore conserved *in expectation*, avoiding systematic bias.
	! ------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, nLocalDepthLayers, nDepthLayersToSeedPhyto, iTimeStep
	real*8, intent(in) :: gridCellArea, NPP
	real*8, dimension(nLocalDepthLayers), intent(in) :: ztop, zbot,  OmegaCalc, TempC, Rho, DynVisco
	real*8, dimension(nDepthLayersToSeedPhyto,nPfts), intent(in) :: probPftProfile
	real*8, dimension(nPfts), intent(in) :: probPftMean
	logical, dimension(nLocalDepthLayers), intent(in) :: validSeedingDepthIdx
	integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms,nDepthLayers), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	integer, dimension(nNewPhytoClusters), intent(out) :: newPhytoClustersIdxs, newPhytoClustersPfts
	real*8, intent(out) :: molesNewOrgCarbon
	
	integer :: iCluster, iPhytoCluster, iRandDepthLayer, iRandPft, iDl, iCompressedDepth, iValid, nCreated
	real*8 :: totalOrgC, availOrgCperCluster, nParticlesPerCluster, randDepth, harvest, orgCperCell, &
		tempScalingFactor, volume, addDiatBiomass, addFlagelBiomass, addCoccoBiomass, &
		addPicoBiomass, remainder, expectedNumParticlesPerCluster, realisedNumParticlesPerCluster
	real*8, dimension(nDepthLayers) :: unusedOrgCByDepth
	logical :: isDiatSelect, isFlagelSelect, isCoccoSelect, isPicoSelect

	! ------ Initialise these (default output) ------
	newPhytoClustersIdxs(:) = 0
    newPhytoClustersPfts(:) = 0
    molesNewOrgCarbon = 0d0
    
    ! ------ Input sanity ------
	totalOrgC = NPP * gridCellArea * timeStep ! mol
	if (.not. ieee_is_finite(totalOrgC) .or. totalOrgC <= 0d0) return
	availOrgCperCluster = totalOrgC / dble(nNewPhytoClusters) ! mol
	if (.not. ieee_is_finite(availOrgCperCluster) .or. availOrgCperCluster <= 0d0) return
    
	! ------ Initialise these ------
	isDiatSelect = .false.
	isFlagelSelect = .false.
	isCoccoSelect = .false.
	isPicoSelect = .false.
	
	addDiatBiomass = 0d0
	addFlagelBiomass = 0d0
	addCoccoBiomass = 0d0
	addPicoBiomass = 0d0
	
	nCreated = 0
    unusedOrgCByDepth(:) = 0d0
    
	tempScalingFactor = MAXVAL(TempC(:))/25d0 ! to scale carbon quota selection (and therefore phytoplankton size) up or down
	
	! ------ Seeding loop ------
	do iPhytoCluster = 1, nNewPhytoClusters

		! ------ Type of phytoplankton functional type (PFT) ------
		iRandPft = RandProbCase(nPfts, probPftMean(:))	
		
		! ------ Depth of the PFT ------	
		select case (choicePhytoSeedingDepthDistrib)	
		case (1) ! SLAMS1.0 - fixed depth distribution
			call random_number(harvest)
			randDepth = zbot(nLocalDepthLayers)*EXP(-5d0*harvest)
			iRandDepthLayer = FindDepthLayerIndex(randDepth, nLocalDepthLayers, ztop(:), zbot(:))
		
		case (2) ! SLAMS2.0 - depth distribution that is environment dependent
			! Sample compressed index from the PDF (1..nDepthLayersToSeedPhyto)
			iCompressedDepth = RandProbCase(nDepthLayersToSeedPhyto, probPftProfile(:,iRandPft))
	
			! Map compressed index -> real depth-layer index (1..nLocalDepthLayers)
			iValid = 0
			iRandDepthLayer = -1
			do iDl = 1, nLocalDepthLayers
				if (.not. validSeedingDepthIdx(iDl)) cycle
				iValid = iValid + 1
				if (iValid == iCompressedDepth) then
					iRandDepthLayer = iDl
					exit
				end if
			end do
			
			if (iRandDepthLayer < 1) then
    			write(*,*) 'ERROR: failed to map compressed depth index to a valid depth layer'
    			call WriteStatusAndStop( )
			end if

			! Obtain the actual depth
        	randDepth = FindDepthLayerMidpointDepth(iRandDepthLayer, nLocalDepthLayers, ztop(:), zbot(:))   
		end select
		
		! ------ Quota per cell ------
		call random_number( harvest )
		if (tempScalingFactor < 0.5d0) then ! cold temperatures, bigger sizes
			harvest = MIN((harvest + 0.20d0), 1d0)
		elseif (tempScalingFactor >= 0.5d0) then ! warm temperatures, smaller sizes
			harvest = MAX((harvest - 0.20d0), 0d0)
		end if
		
		orgCperCell = C_quota_pft_min(iRandPft) + harvest*(C_quota_pft_max(iRandPft)-C_quota_pft_min(iRandPft)) 
		if (availOrgCperCluster < orgCperCell) then ! this shouldn't happen after all prior checks, but just in case
			unusedOrgCByDepth(iRandDepthLayer) = unusedOrgCByDepth(iRandDepthLayer) + availOrgCperCluster
			cycle 
		end if

		! ------ No. particles per cluster ------
		! nParticlesPerCluster is real*8 but always integer-valued
		! Impose discreteness on phytoplankton particles using stochastic rounding.
		! The expected number of particles per cluster is generally non-integer because
		! available organic carbon is continuous. To avoid systematic carbon loss (floor)
		! or gain (round), apply stochastic rounding
		expectedNumParticlesPerCluster = SafeDivide(availOrgCperCluster, orgCperCell, 0d0)
		if (expectedNumParticlesPerCluster <= 0d0) then
			unusedOrgCByDepth(iRandDepthLayer) = unusedOrgCByDepth(iRandDepthLayer) + availOrgCperCluster
			cycle
		end if

		realisedNumParticlesPerCluster = SafeFloorNonNegative(expectedNumParticlesPerCluster) ! round down
		remainder = expectedNumParticlesPerCluster - realisedNumParticlesPerCluster ! used as the probability of adding one extra particle
		remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
		
		call random_number(harvest)
		if (harvest < remainder) then
			nParticlesPerCluster = realisedNumParticlesPerCluster + 1d0
		else
			nParticlesPerCluster = realisedNumParticlesPerCluster
		end if
		if (nParticlesPerCluster < 1d0) then
    		write(*,*)  'ERROR: nParticlesPerCluster < 1 in PhytoplanktonProduction'
    		call WriteStatusAndStop()
		end if
		
		! ------ Assign next free slot ------
		if (iLastLocus + 1 > nClusters) then
            write(*,*) 'ERROR: particle array exhausted in PhytoplanktonProduction'
            call WriteStatusAndStop()
        end if
		iCluster = iLastLocus + 1

		! ------ Compute the attributes ------ 	 	
		call ComputeParticleAttributes(particle, nClusters, iCluster, 3, iRandPft, 1, 0, & 
			iTimeStep, OmegaCalc(iRandDepthLayer), randDepth, nParticlesPerCluster, &
			orgCperCell, Rho(iRandDepthLayer), DynVisco(iRandDepthLayer)) ! -- in calcparticleattributes.F90
		
		! ------ Store diagnostics ------
		SMSterm(iPrimProdOrgC,iRandDepthLayer) = SMSterm(iPrimProdOrgC,iRandDepthLayer) &
			+ particle(iCluster)%molesOrgC*nParticlesPerCluster ! mol orgC
		SMSterm(iPrimProdCaCO3,iRandDepthLayer) = SMSterm(iPrimProdCaCO3,iRandDepthLayer) &
			+ particle(iCluster)%molesMineral(iCalcite)*nParticlesPerCluster ! mol CaCO3
		SMSterm(iPrimProdOpal,iRandDepthLayer) = SMSterm(iPrimProdOpal,iRandDepthLayer) &
			+ particle(iCluster)%molesMineral(iOpal)*nParticlesPerCluster ! mol opal
	
		if (iRandPft == 1) then
			auxTerm(iFreshDiatCellQuota) = auxTerm(iFreshDiatCellQuota) + orgCperCell
			auxCount(iFreshDiatCellQuota) = auxCount(iFreshDiatCellQuota) + 1
			addDiatBiomass = addDiatBiomass + particle(iCluster)%molesOrgC*nParticlesPerCluster ! mol C
			isDiatSelect = .true.
		elseif (iRandPft == 2) then
			auxTerm(iFreshFlagelCellQuota) = auxTerm(iFreshFlagelCellQuota) + orgCperCell
			auxCount(iFreshFlagelCellQuota) = auxCount(iFreshFlagelCellQuota) + 1
			addFlagelBiomass = addFlagelBiomass + particle(iCluster)%molesOrgC*nParticlesPerCluster ! mol C
			isFlagelSelect = .true.
		elseif (iRandPft == 3) then			
			auxTerm(iFreshCoccoCellQuota) = auxTerm(iFreshCoccoCellQuota) + orgCperCell
			auxCount(iFreshCoccoCellQuota) = auxCount(iFreshCoccoCellQuota) + 1
			addCoccoBiomass = addCoccoBiomass + particle(iCluster)%molesOrgC*nParticlesPerCluster ! mol C
			isCoccoSelect = .true.
		elseif (iRandPft == 4) then	
			auxTerm(iFreshPicoCellQuota) = auxTerm(iFreshPicoCellQuota) + orgCperCell
			auxCount(iFreshPicoCellQuota) = auxCount(iFreshPicoCellQuota) + 1
			addPicoBiomass = addPicoBiomass + particle(iCluster)%molesOrgC*nParticlesPerCluster ! mol C
			isPicoSelect = .true.
		end if
		
 		! ------ Update on success seeding ------
		iLastLocus = iLastLocus + 1
		nCreated = nCreated + 1
		newPhytoClustersIdxs(nCreated) = iCluster
        newPhytoClustersPfts(nCreated) = iRandPft
        molesNewOrgCarbon = molesNewOrgCarbon + particle(iCluster)%molesOrgC*nParticlesPerCluster	
		
	end do
	
	! ------ Redirect to microbial pool, so mass is not lost ------
	if (SUM(unusedOrgCByDepth) > 0d0) then
		do iDl = 1, nDepthLayers
			if (unusedOrgCByDepth(iDl) > 0d0) then
				SMSterm(iMicrobSolubOrgC,iDl) = SMSterm(iMicrobSolubOrgC,iDl) + unusedOrgCByDepth(iDl)
			end if
		end do
	end if
	
	! ------ Store diagnostics ------
	volume = 0d0
	do iDl = 1, nDepthLayersToSeedPhyto
		if (.not. validSeedingDepthIdx(iDl)) cycle
		volume = volume + (zbot(iDl)-ztop(iDl))*gridCellArea ! 1 m3
	end do

	if (isDiatSelect) then
		auxTerm(iDiatBiomass) = auxTerm(iDiatBiomass) + addDiatBiomass/volume ! mol C m-3
		auxCount(iDiatBiomass) = auxCount(iDiatBiomass) + 1
	end if
	if (isFlagelSelect) then	
		auxTerm(iFlagelBiomass) = auxTerm(iFlagelBiomass) + addFlagelBiomass/volume ! mol C m-3
		auxCount(iFlagelBiomass) = auxCount(iFlagelBiomass) + 1	
	end if
	if (isCoccoSelect) then
		auxTerm(iCoccoBiomass) = auxTerm(iCoccoBiomass) + addCoccoBiomass/volume ! mol C m-3
		auxCount(iCoccoBiomass) = auxCount(iCoccoBiomass) + 1
	end if
	if (isPicoSelect) then
		auxTerm(iPicoBiomass) = auxTerm(iPicoBiomass) + addPicoBiomass/volume ! mol C m-3
		auxCount(iPicoBiomass) = auxCount(iPicoBiomass) + 1
	end if

end subroutine PhytoplanktonProduction

! ========================================================================================

subroutine TepProduction(particle, nClusters, iLastLocus, nLocalDepthLayers, ztop, zbot, &
	SMSterm, newPhytoClustersIdxs, newPhytoClustersPfts, Rho, DynVisco, iTimeStep)

	! ------------------------------------------------------------------------------------
	! TEP are produced out of the organic matter of large PFTs (diatoms and flagellated 
	! phytoplankton)
	! ------------------------------------------------------------------------------------	

	integer, intent(in) :: nClusters, nLocalDepthLayers, iTimeStep
	integer, dimension(nNewPhytoClusters), intent(in) :: newPhytoClustersIdxs, newPhytoClustersPfts
	real*8, dimension(nLocalDepthLayers), intent(in) :: ztop, zbot, Rho, DynVisco
	integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms,nDepthLayers), intent(inout) :: SMSterm
	
  	integer :: nSelected, nLarge, iPhytoCluster, iDl, iTepCluster, iCluster
	integer, allocatable :: selectedPhyto(:)
  	real*8 :: releasedOrgC, depthPhyto, randTepCquota, nTepParticlesPerCluster, &
  		availOrgCperCluster, expectedNumTepParticles, realisedNumTepParticles, remainder, &
  		harvest
  	real*8, dimension(nDepthLayers) :: unusedOrgCByDepth
  	real*8, allocatable :: u(:)
  	
	! ------ Select diatoms & flagellates once ------
  	nLarge = count( (newPhytoClustersPfts == 1) .or. (newPhytoClustersPfts == 2) )
  	if (nLarge == 0) return

	allocate(selectedPhyto(nLarge))
  	nSelected = 0
	releasedOrgC = 0d0

	! ------ Exudation + immediate shrink/recompute for the producing phyto ------
	do iPhytoCluster = 1, nNewPhytoClusters
		if (newPhytoClustersPfts(iPhytoCluster) == 1 .or. newPhytoClustersPfts(iPhytoCluster) == 2) then
      		nSelected = nSelected + 1
      		selectedPhyto(nSelected) = newPhytoClustersIdxs(iPhytoCluster)

			releasedOrgC = releasedOrgC + phyto_exudation_frac * particle(selectedPhyto(nSelected))%molesOrgC &
            	* particle(selectedPhyto(nSelected))%nPxC

			! Shrink orgC in cell
      		particle(selectedPhyto(nSelected))%molesOrgC = &
    			particle(selectedPhyto(nSelected))%molesOrgC * (1d0 - phyto_exudation_frac)

			depthPhyto  = particle(selectedPhyto(nSelected))%depth
      		iDl = FindDepthLayerIndex(depthPhyto, nLocalDepthLayers, ztop, zbot)
			
			! Recompute mass/geom/density/velocity (porosity & stickiness unchanged)				
			call ParticleDryMass(particle, nClusters, selectedPhyto(nSelected))					
			call ParticleMaterialVolume(particle, nClusters, selectedPhyto(nSelected))
			call ParticleRadius(particle, nClusters, selectedPhyto(nSelected))
			call ParticleDensity(particle, nClusters, selectedPhyto(nSelected), Rho(iDl)) 
			call ParticleSettlingVelocity(particle, nClusters, selectedPhyto(nSelected), Rho(iDl), DynVisco(iDl) ) ! m d-1	
			call CheckParticleSanity(particle, nClusters, selectedPhyto(nSelected), 'phytoplankton cell after exudation')				
		end if
		
	end do
	
	! ------ Safety checks ------
	if (.not. ieee_is_finite(releasedOrgC) .or. releasedOrgC <= 0d0) then
    	deallocate(selectedPhyto)
    	return
  	end if
  	allocate(u(nNewPhytoTepClusters)); call random_number(u)
  	availOrgCperCluster = releasedOrgC / dble(nNewPhytoTepClusters)
  	if (.not. ieee_is_finite(availOrgCperCluster) .or. availOrgCperCluster <= 0d0) then
    	deallocate(u, selectedPhyto)
    	return
	end if
	
	! Initialise this
  	unusedOrgCByDepth(:) = 0d0
  	
  	! ------ Create TEPs ------
	do iTepCluster = 1, nNewPhytoTepClusters	

		! ------ Depth: reuse a producer if enough, else fallback to exponential ------
		if (iTepCluster <= nSelected) then
      		depthPhyto = particle(selectedPhyto(iTepCluster))%depth
    	else
      		depthPhyto = zbot(nLocalDepthLayers)*exp(-5d0*u(iTepCluster))
    	end if
    	iDl = FindDepthLayerIndex(depthPhyto, nLocalDepthLayers, ztop, zbot)				

		! ------ Quota per cell ------
    	randTepCquota = C_quota_TEP_min + (C_quota_TEP_max - C_quota_TEP_min)*u(iTepCluster)
    	if (availOrgCperCluster < randTepCquota) then ! this shouldn't happen, but just in case
    		unusedOrgCByDepth(iDl) = unusedOrgCByDepth(iDl) + availOrgCperCluster
			cycle 
		end if
			
 		! ------ Compute nTepParticlesPerCluster, a real*8 but integer-valued ------
 		expectedNumTepParticles = SafeDivide(availOrgCperCluster, randTepCquota, 0d0)
		if (expectedNumTepParticles <= 0d0) then
			unusedOrgCByDepth(iDl) = unusedOrgCByDepth(iDl) + availOrgCperCluster
			cycle
		end if
		realisedNumTepParticles = SafeFloorNonNegative(expectedNumTepParticles) ! round down
		remainder = expectedNumTepParticles - realisedNumTepParticles ! used as the probability of adding one extra particle
		remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
		
		call random_number(harvest)
		if (harvest < remainder) then
			nTepParticlesPerCluster = realisedNumTepParticles + 1d0
		else
			nTepParticlesPerCluster = realisedNumTepParticles
		end if
		if (nTepParticlesPerCluster < 1d0) then
    		write(*,*) 'ERROR: nTepParticlesPerCluster < 1 in TepProduction'
    		call WriteStatusAndStop()
		end if
		
		! ------ Assign next free slot ------
		if (iLastLocus + 1 > nClusters) then
            write(*,*) 'ERROR: particle array exhausted in TepProduction'
            call WriteStatusAndStop()
        end if
		iCluster = iLastLocus + 1
		
		! ------ Compute the attributes ------ 			
		call ComputeParticleAttributes(particle, nClusters, iCluster, 5, 0, 0, 0, iTimeStep, &
			0d0, depthPhyto, nTepParticlesPerCluster, randTepCquota, Rho(iDl), DynVisco(iDl))		

		! ------ Store diagnostics ------
		SMSterm(iProdTepPhyto,iDl) = SMSterm(iProdTepPhyto,iDl) + particle(iCluster)%molesTepC*nTepParticlesPerCluster ! mol TEP
			
		! ------ Update on success seeding ------
		iLastLocus = iLastLocus + 1
	end do
	
	! ------ Redirect to microbial pool, so mass is not lost ------
	if (SUM(unusedOrgCByDepth) > 0d0) then
		do iDl = 1, nDepthLayers
			if (unusedOrgCByDepth(iDl) > 0d0) then
				SMSterm(iMicrobSolubTepC,iDl) = SMSterm(iMicrobSolubTepC,iDl) + unusedOrgCByDepth(iDl)
			end if
		end do
	end if

	deallocate(u,selectedPhyto)

end subroutine TepProduction

! ========================================================================================

subroutine ClayDeposition(particle, nClusters, iLastLocus, gridCellArea, SMSterm, &
    accumAeolClay, waterRhoSurface, dynViscoSurface, iTimeStep)

	integer, intent(in) :: nClusters, iTimeStep
	real*8, intent(in) :: gridCellArea, accumAeolClay, waterRhoSurface, dynViscoSurface
	integer, intent(inout) :: iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
  	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
 	
 	integer :: iClayCluster, iCluster
  	real*8 :: availClayPerCluster, expectedNumClayParticles, realisedNumClayParticles, remainder, &
  		nClayParticlesPerCluster, clayDepth, harvest

  	availClayPerCluster = accumAeolClay*gridCellArea/(MOLAR_MASS_CLAY*dble(nNewClayClusters)) ! mol
  	if (.not. ieee_is_finite(availClayPerCluster) .or. availClayPerCluster <= 0d0) return
  	
  	! If even one clay particle cannot be formed per intended cluster, redirect ALL mass once
  	if (availClayPerCluster < clay_quota) then
  		SMSterm(iMicrobSolubClay) = SMSterm(iMicrobSolubClay) + accumAeolClay*gridCellArea/MOLAR_MASS_CLAY ! fallback --> mass is not lost, it is redirected
  		return
  	end if
  	
  	! ------ Compute nClayParticlesPerCluster, a real*8 but integer-valued ------
	expectedNumClayParticles = SafeDivide(availClayPerCluster, clay_quota, 0d0)
	if (expectedNumClayParticles <= 0d0) return ! shouldn't happen after all the checks before, but just in case...
	realisedNumClayParticles = SafeFloorNonNegative(expectedNumClayParticles) ! round down
	remainder = expectedNumClayParticles - realisedNumClayParticles ! used as the probability of adding one extra particle
	remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety

	do iClayCluster = 1, nNewClayClusters
		
		call random_number(harvest)
		if (harvest < remainder) then
			nClayParticlesPerCluster = realisedNumClayParticles + 1d0
		else
			nClayParticlesPerCluster = realisedNumClayParticles
		end if
		if (nClayParticlesPerCluster < 1d0) then
    		write(*,*) 'INTERNAL ERROR: nClayParticlesPerCluster < 1 in ClayDeposition'
   	 		call WriteStatusAndStop()
		end if
		
		! ------ Assign next free slot ------
		if (iLastLocus + 1 > nClusters) then
            write(*,*) 'ERROR: particle array exhausted in ClayDeposition'
            call WriteStatusAndStop()
        end if
		iCluster = iLastLocus + 1
		
		! ------ Compute attributes ------
		clayDepth = 0d0 ! set to very top surface
		call ComputeParticleAttributes(particle, nClusters, iCluster, 6, 0, 0, 0, iTimeStep, &
			0d0, clayDepth, nClayParticlesPerCluster, 0d0, waterRhoSurface, dynViscoSurface)

		SMSterm(iDepoClay) = SMSterm(iDepoClay) + particle(iCluster)%molesMineral(iClay)*nClayParticlesPerCluster
		
		! ------ Update on success seeding ------
		iLastLocus = iLastLocus + 1
	end do

end subroutine ClayDeposition

! ========================================================================================

subroutine FixedProduction(particle, nClusters, iLastLocus, SMSterm, iTimeStep)

	integer, intent(in) :: nClusters, iTimeStep
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle	
	integer, intent(inout) :: iLastLocus
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	
	integer :: iFixedCluster, iCluster, firstNewCluster, lastNewCluster
	integer, dimension(nNewFixedClusters) :: newFixedClusters
	real*8 :: nParticlesPerCluster, orgCperCell, clusterDepth, waterRho, waterKinVisco, waterDynVisco
	integer :: i

	firstNewCluster = iLastLocus + 1
	lastNewCluster = iLastLocus + nNewFixedClusters
	newFixedClusters(:) = [(i, i = firstNewCluster, lastNewCluster, 1)]

	do iFixedCluster = 1, nNewFixedClusters
		iCluster = newFixedClusters(iFixedCluster)
		clusterDepth = 0.10d0 ! m
		waterRho = 1.0275d0 ! g cm-3	
		waterKinVisco = 1d-6 ! m2 s-1
		waterDynVisco = 1d4*waterKinVisco * waterRho ! g cm-1 s-1
		nParticlesPerCluster = 1d8 ! (1d8 x 1d4 particles) / 1d4 clusters (1d12 particles in 1d6 cm3 is 1d6 particles/cm3, as Tinna's and Jackson)
		orgCperCell = 1.67d-4 * (20d-4/2d0)**2.28d0 ! mol C per cell
			
		call ComputeParticleAttributes(particle, nClusters, iCluster, 3, 4, 1, 0, iTimeStep, &
			0d0, clusterDepth, nParticlesPerCluster, orgCperCell, waterRho, waterDynVisco) ! create picophytoplankton cells
		
		SMSterm(iPrimProdOrgC) = SMSterm(iPrimProdOrgC) + particle(iCluster)%molesOrgC*nParticlesPerCluster ! mol orgC

		if (iFixedCluster == 1) then
			write(*,*) 'Init type attributes'
			write(*,*) '  diameter:', 2d0*particle(iCluster)%radius
			write(*,*) '  volume  :', particle(iCluster)%solidVolume
			write(*,*) '  velocity:', particle(iCluster)%velocity
			write(*,*) '  density :', particle(iCluster)%density
			write(*,*) '  nPxC    :', particle(iCluster)%nPxC
			write(*,*) '  C, mole :', particle(iCluster)%molesOrgC
		end if

	end do

	iLastLocus = lastNewCluster
	
end subroutine FixedProduction

! ========================================================================================

end module injectsurfaceparticles