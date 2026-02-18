#include "blockdefinitions.h"

module modeloutput

! ----------------------------------------------------------------------------------------
! This module handles writing output files as well as printing information to the log file.
! ----------------------------------------------------------------------------------------

use particlestructure, only: lagrangianStateVars
use modelcounters, only: particleCollectionPeriod, iLastSliceSedTrap, iLastSliceAvgAtt, &
	iLastSliceAvgVolAtt, iLastSliceEulerian, iLastSliceAux, iLastSliceStats, iLastSliceSnapshot, &
	iLastSliceParticleNumSnapshot, nCollisionsType1, nCollisionsType2, nCollisionsType3, &
	nCollisionsType4, nTimesEnteringCollLoop, nZooIngestionEvents, nZooFragmentationEvents, &
	nParticlesEvaluatedForGrazing, nFaecalPelletClustersProduced, nZooDeadClustersProduced, &
	nMicrobRespiredClusters, nFragmentedBigClusters, nRespiredTinyClusters, nMineralClustersDissolved, &
	nClustersPhotolysed, nClustersAtSeafloor
use modelconstants, only: MOLAR_MASS_CARBON, MOLAR_MASS_CACO3, MOLAR_MASS_OPAL, &
	MOLAR_MASS_CLAY, SECONDS_PER_DAY
use modelparameters, only: maxNumSmsTerms, maxNumAuxTerms, nDepthLayers, nSizeClasses, &
	nVeloClasses, nVolumeClasses, nMainParticleTypes, particleSizeClasses, particleVeloClasses, &
	particleVolumeClasses, iOpal, iCalcite, iClay, iZooRespOrgC, iZooRespTepC, &
	iZooIngestOrgC, iZooIngestTepC, iZooEgestOrgC, iZooEgestTepC, iZooExcretOrgC, &
	iZooExcretTepC, iZooSolubOrgC, iZooSolubTepC, iZooDissolCaCO3, iPrimProdOrgC, &
	iProdTepPhyto, iPrimProdCaCO3, iPrimProdOpal, iZooDeathOrgC, iDepoClay, iMicrobRespOrgC, &
	iMicrobRespTepC, iMicrobSolubOrgC, iMicrobSolubTepC, iZooNumber, iZooBiomass, &
	iNightDvmUpperBound, iNightDvmLowerBound, iDayDvmUpperBound, iDayDvmLowerBound, iDissolCaCO3, &
	iDissolOpal, iFreshDiatCellQuota, iFreshFlagelCellQuota, iFreshCoccoCellQuota, iFreshPicoCellQuota, &
	iZooSpecRespRate, iMicrobSpecRespRate, iCoagulationSuccess, iDiatBiomass, iFlagelBiomass, &
	iCoccoBiomass, iPicoBiomass, iEuphoticDepth, iCollisionKernel, iNumParticlesEvalEncounter, &
	iNumParticlesEncountered, iAvgProbDiat, iAvgProbFlagel, iAvgProbCocco, iAvgProbPico, &
	iBrownianKernel, iShearKernel, iSettlingKernel, detection_limit_poc, &
	filenameFlux, filenameFluxSf, filenameAvgAttSizeClass, filenameAvgAttVeloClass, &
	filenameAvgAttVolumeClass, filenameAvgAttSizeClassSf, filenameAvgAttVeloClassSf, &
	filenameInstAvgAttSizeClass, filenameInstAvgAttVeloClass, filenameInstAvgAttSizeClassSf, &
	filenameInstAvgAttVeloClassSf, filenameAvgAttMainType, filenameAvgAttMainTypeSf, filenameInstAvgAttMainType, &
	filenameInstAvgAttMainTypeSf, filenameLossTerms, filenameSms, filenameSmsIntegrated, filenameAuxTerms, &
	filenameStatsNumClusters, filenameStatsNumParticles, filenameControl, filenameClustersInteger, &
	filenameClustersReal
use modelgrid, only: gridCellArea, ztop, zbot
use modeleulerianvariables, only: accumSMS, avgSMSflux, avgSMSfluxIntegrated, auxTerm, auxCount
use modelparticlecollection, only: nSedTrapDeployDepths, nImagingDeployDepths, sedTrap, sedTrapSf, &
	cameraUpperBoundaries, cameraLowerBoundaries, attsInSizeClass, attsInVeloClass, attsInVolumeClass, &
	attsInSizeClassSf, attsInVeloClassSf, attsInSizeClassInst, attsInVeloClassInst, &
	attsInSizeClassSfInst, attsInVeloClassSfInst, particleFlux, particleFluxSf, &
	avgAttsInSizeClass, avgAttsInVeloClass, avgAttsInVolumeClass, avgAttsInSizeClassSf, &
	avgAttsInVeloClassSf, avgAttsInSizeClassInst, avgAttsInVeloClassInst, avgAttsInSizeClassSfInst, &
	avgAttsInVeloClassSfInst, attsInMainType, attsInMainTypeSf, attsInMainTypeInst, attsInMainTypeSfInst, &
	avgAttsInMainType, avgAttsInMainTypeSf, avgAttsInMainTypeInst, avgAttsInMainTypeSfInst, &
	cameraSampledVolume, lossTerms
use sanitychecks, only: PrintParticleProperties, WriteStatusAndStop
use findfunctions, only: FindParticleAttributeClass, FindParticleVolumeClass, FindDepthLayerIndex
	
implicit none
external :: write_r8_field, write_i4_field, write_i8_field
private
public :: GetParticleAverageAttributes, WriteParticleFlux, WriteParticleAverageAttributes, &
	WriteEulerianVariables, WriteAuxiliaryTerms, WriteAndPrintModelStatistics, &
	WriteSnapshotsOfParticleAttributes, WriteModelClosureInformation, WriteInstantaneousSnapshots, &
	GetParticleAverageAttributesByVolumeClass, WriteTestSmoluchowski
integer :: iProfileMonitor = 1

contains

! ========================================================================================

subroutine GetParticleAverageAttributes(nZ, sizeClass, veloClass, mainType, particle, nClusters, &
	iCluster, iYear)

	!-------------------------------------------------------------------------------------
	! The philosophy of this subroutine is to visualise how standard particle size categories
	! (and standard velocity categories) look in terms of density, sinking speed (size), porosity, 
	! stickiness and type of material content. It fills up the arrays sizeClass and veloClass 
	! by accumulating particle attribute values. These arrays are then used in the subroutine 
	! WriteParticleAverageAttributes (in this module, further down) where the average is computed.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nZ, nClusters, iCluster, iYear	
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	real*8, dimension(nSizeClasses,nZ,14), intent(inout) :: sizeClass ! can be the water column camera (nZ=nLocalDepthLayers) or the seafloor camera (nZ=1)
	real*8, dimension(nVeloClasses,nZ,14), intent(inout) :: veloClass ! can be the water column camera (nZ=nLocalDepthLayers) or the seafloor camera (nZ=1)
	real*8, dimension(nMainParticleTypes,nZ,7), intent(inout) :: mainType ! can be the water column camera (nZ=nLocalDepthLayers) or the seafloor camera (nZ=1)
			
	integer :: iSc, iVc, iDh, iTy
	real*8 :: nParticlesPerCluster

	iSc = 0
	iVc = 0
	iDh = 0
	iTy = 0

	! Find the size class, velo class, depth horizon and overall type of the particle
	iSc = FindParticleAttributeClass( nSizeClasses, particleSizeClasses(:), particle(iCluster)%radius*2d0 )
	iVc = FindParticleAttributeClass( nVeloClasses, particleVeloClasses(:), particle(iCluster)%velocity )
	if (particle(iCluster)%phase == 1) then ! particle in the water column
		iDh = FindDepthLayerIndex(particle(iCluster)%depth, nImagingDeployDepths, &
			cameraUpperBoundaries(:), cameraLowerBoundaries(:))
	else if (particle(iCluster)%phase == 2) then ! particle on the seafloor
		iDh = 1
	end if
	iTy = particle(iCluster)%initType

	if (iDh == 0 .or. iSc == 0 .or. iVc == 0 .or. iTy == 0) then
		call PrintParticleProperties(particle, nClusters, iCluster, 'ERROR: in GetParticleAverageAttributes, something is wrong')
    	call WriteStatusAndStop( )
	end if
	
	nParticlesPerCluster = particle(iCluster)%nPxC
	mainType(iTy,iDh,1) = mainType(iTy,iDh,1) + nParticlesPerCluster
	mainType(iTy,iDh,2) = mainType(iTy,iDh,2) + particle(iCluster)%mass*nParticlesPerCluster
	mainType(iTy,iDh,3) = mainType(iTy,iDh,3) + (particle(iCluster)%molesOrgC+particle(iCluster)%molesTepC)*nParticlesPerCluster				
	mainType(iTy,iDh,4) = mainType(iTy,iDh,4) + 2d0*particle(iCluster)%radius*nParticlesPerCluster
	mainType(iTy,iDh,5) = mainType(iTy,iDh,5) + particle(iCluster)%density*nParticlesPerCluster
	mainType(iTy,iDh,6) = mainType(iTy,iDh,6) + particle(iCluster)%velocity*nParticlesPerCluster
	mainType(iTy,iDh,7) = mainType(iTy,iDh,7) + particle(iCluster)%porosity*nParticlesPerCluster
	
	! Only record organics that are sinking
	if ((particle(iCluster)%molesOrgC + particle(iCluster)%molesTepC) > 0d0 &
		.and. particle(iCluster)%velocity > 0d0) then

		! Num of particles per cluster - 1
		sizeClass(iSc,iDh,1) = sizeClass(iSc,iDh,1) + nParticlesPerCluster
		veloClass(iVc,iDh,1) = veloClass(iVc,iDh,1) + nParticlesPerCluster

		! Density - 2
		sizeClass(iSc,iDh,2) = sizeClass(iSc,iDh,2) + particle(iCluster)%density*nParticlesPerCluster
		veloClass(iVc,iDh,2) = veloClass(iVc,iDh,2) + particle(iCluster)%density*nParticlesPerCluster

		! Velocity/radius - 3
		! Attention: record velocity for the size classes, record size for the velocity classes	
		sizeClass(iSc,iDh,3) = sizeClass(iSc,iDh,3) + particle(iCluster)%velocity*nParticlesPerCluster
		veloClass(iVc,iDh,3) = veloClass(iVc,iDh,3) + 2d0*particle(iCluster)%radius*nParticlesPerCluster

		! Stickiness - 4 			
		sizeClass(iSc,iDh,4) = sizeClass(iSc,iDh,4) + particle(iCluster)%stickiness*nParticlesPerCluster
		veloClass(iVc,iDh,4) = veloClass(iVc,iDh,4) + particle(iCluster)%stickiness*nParticlesPerCluster

		! Porosity - 5		
		sizeClass(iSc,iDh,5) = sizeClass(iSc,iDh,5) + particle(iCluster)%porosity*nParticlesPerCluster		
		veloClass(iVc,iDh,5) = veloClass(iVc,iDh,5) + particle(iCluster)%porosity*nParticlesPerCluster

		! Solid volume - 6
		sizeClass(iSc,iDh,6) = sizeClass(iSc,iDh,6) + particle(iCluster)%solidVolume*nParticlesPerCluster
		veloClass(iVc,iDh,6) = veloClass(iVc,iDh,6) + particle(iCluster)%solidVolume*nParticlesPerCluster
	
		! Fractal dimension - 7
		sizeClass(iSc,iDh,7) = sizeClass(iSc,iDh,7) + particle(iCluster)%fracDim*nParticlesPerCluster
		veloClass(iVc,iDh,7) = veloClass(iVc,iDh,7) + particle(iCluster)%fracDim*nParticlesPerCluster

		! Radius PP - 8
		sizeClass(iSc,iDh,8) = sizeClass(iSc,iDh,8) + 2d0*particle(iCluster)%radiusPp*nParticlesPerCluster			
		veloClass(iVc,iDh,8) = veloClass(iVc,iDh,8) + 2d0*particle(iCluster)%radiusPp*nParticlesPerCluster

		! Depth - 9
		sizeClass(iSc,iDh,9) = sizeClass(iSc,iDh,9) + particle(iCluster)%depth*nParticlesPerCluster			
		veloClass(iVc,iDh,9) = veloClass(iVc,iDh,9) + particle(iCluster)%depth*nParticlesPerCluster		

		! Moles orgC - 10			
		sizeClass(iSc,iDh,10) = sizeClass(iSc,iDh,10) + particle(iCluster)%molesOrgC*nParticlesPerCluster
		veloClass(iVc,iDh,10) = veloClass(iVc,iDh,10) + particle(iCluster)%molesOrgC*nParticlesPerCluster	

		! Moles tepC - 11	
		sizeClass(iSc,iDh,11) = sizeClass(iSc,iDh,11) + particle(iCluster)%molesTepC*nParticlesPerCluster
		veloClass(iVc,iDh,11) = veloClass(iVc,iDh,11) + particle(iCluster)%molesTepC*nParticlesPerCluster

		! Moles opal - 12				
		sizeClass(iSc,iDh,12) = sizeClass(iSc,iDh,12) + particle(iCluster)%molesMineral(iOpal)*nParticlesPerCluster		
		veloClass(iVc,iDh,12) = veloClass(iVc,iDh,12) + particle(iCluster)%molesMineral(iOpal)*nParticlesPerCluster		

		! Moles calcite - 13			
		sizeClass(iSc,iDh,13) = sizeClass(iSc,iDh,13) + particle(iCluster)%molesMineral(iCalcite)*nParticlesPerCluster	
		veloClass(iVc,iDh,13) = veloClass(iVc,iDh,13) + particle(iCluster)%molesMineral(iCalcite)*nParticlesPerCluster			
										  
		! Moles clay - 14
		sizeClass(iSc,iDh,14) = sizeClass(iSc,iDh,14) + particle(iCluster)%molesMineral(iClay)*nParticlesPerCluster
		veloClass(iVc,iDh,14) = veloClass(iVc,iDh,14) + particle(iCluster)%molesMineral(iClay)*nParticlesPerCluster
			
	end if
						
end subroutine GetParticleAverageAttributes

! ========================================================================================

subroutine WriteParticleAverageAttributes(nAveragingTimeSteps,localSeafloorDepth,nProfiles)

	!-------------------------------------------------------------------------------------
	! Calculate the average attribute for that particle size class or velo class by dividing 
	! the accumulated attribute by the number of accumulated particles. This subroutine has 
	! two parts: first average particle attributes are calculated, second the particle numebr
	! slice is populated.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nProfiles, nAveragingTimeSteps
	real*8, intent(in) :: localSeafloorDepth
	
	integer :: iDh, iSc, iVc, iAtt, iTy, newAvgAttsSlice, iProfile, i
	real*8 :: nParticles20m, nParticles50m, nParticles100m, nParticles500m, &
		nParticles1km, mass20m, mass50m, mass100m, mass500m, mass1km, sizeClassRange
    character(len=30) :: fmt1, fmt2
	
	avgAttsInSizeClass(:,:,:,:) 	= 0d0
	avgAttsInVeloClass(:,:,:,:) 	= 0d0
	avgAttsInSizeClassSf(:,:,:,:)	= 0d0
	avgAttsInVeloClassSf(:,:,:,:) 	= 0d0
	avgAttsInMainType(:,:,:,:) 		= 0d0
	avgAttsInMainTypeSf(:,:,:,:) 	= 0d0
	
	! First, record average particle attributes
	
    do iProfile = 1, nProfiles
    
    	!---------------------------------------------------------------------------------
		! Water column recording
		!---------------------------------------------------------------------------------
		
		do iDh = 1, nImagingDeployDepths
		
			! Particle type classifier
			do iTy = 1, nMainParticleTypes
				avgAttsInMainType(iTy,iDh,1,iProfile) = attsInMainType(iTy,iDh,1,iProfile) ! number of particles
				avgAttsInMainType(iTy,iDh,2,iProfile) = attsInMainType(iTy,iDh,2,iProfile) ! total mass
				avgAttsInMainType(iTy,iDh,3,iProfile) = attsInMainType(iTy,iDh,3,iProfile) ! total POC (orgC+TEPC)
				do iAtt = 4, 7 ! avg. diameter of a particle, avg. density of a particle, avg. velocity of a particle, avg. porosity of a particle
					if (attsInMainType(iTy,iDh,1,iProfile) > 0) then
						avgAttsInMainType(iTy,iDh,iAtt,iProfile) = attsInMainType(iTy,iDh,iAtt,iProfile) &
							/attsInMainType(iTy,iDh,1,iProfile) ! 1=nParticles
					end if				
				end do
			end do		

			do iAtt = 2, 14

		 		! Size classifier in the water column	
				do iSc = 1, nSizeClasses
				  	if (attsInSizeClass(iSc,iDh,1,iProfile) > 0) then
						avgAttsInSizeClass(iSc,iDh,iAtt,iProfile) = attsInSizeClass(iSc,iDh,iAtt,iProfile) &
														 / attsInSizeClass(iSc,iDh,1,iProfile) ! 1=nParticles
				  	end if				
			  	end do
								
		  		! Velocity classifier	in the water column	
			  	do iVc = 1, nVeloClasses	
				  	if (attsInVeloClass(iVc,iDh,1,iProfile) > 0) then 	
					  	avgAttsInVeloClass(iVc,iDh,iAtt,iProfile) = attsInVeloClass(iVc,iDh,iAtt,iProfile) &
														 / attsInVeloClass(iVc,iDh,1,iProfile) ! 1=nParticles
				  	end if	
			  	end do

		  	end do ! iAtt	
				
	  	end do ! iDh

		!---------------------------------------------------------------------------------
		! Seafloor recording
		!---------------------------------------------------------------------------------
	
		! Particle type classifier	
		do iTy = 1, nMainParticleTypes
			do iAtt = 4, 7 ! avg. diameter of a particle, avg. density of a particle, avg. velocity of a particle, avg. porosity of a particle
				if (attsInMainTypeSf(iTy,1,1,iProfile) > 0) then
					avgAttsInMainTypeSf(iTy,1,iAtt,iProfile) = attsInMainTypeSf(iTy,1,iAtt,iProfile) &
						/attsInMainTypeSf(iTy,1,1,iProfile) ! 1=nParticles
				end if				
			end do
			avgAttsInMainTypeSf(iTy,1,1,iProfile) = attsInMainTypeSf(iTy,1,1,iProfile) ! number of particles
			avgAttsInMainTypeSf(iTy,1,2,iProfile) = attsInMainTypeSf(iTy,1,2,iProfile) ! total mass
			avgAttsInMainTypeSf(iTy,1,3,iProfile) = attsInMainTypeSf(iTy,1,3,iProfile) ! total organic carbon
		end do
			
		do iAtt = 2, 14
		
			! Size classifier at seafloor		
			do iSc = 1, nSizeClasses
				if (attsInSizeClassSf(iSc,1,1,iProfile) > 0) then
	  				avgAttsInSizeClassSf(iSc,1,iAtt,iProfile) = attsInSizeClassSf(iSc,1,iAtt,iProfile) &
														 / attsInSizeClassSf(iSc,1,1,iProfile) ! 1=nParticles
				end if
			end do
				
			! Velocity classifier at seafloor			
			do iVc = 1, nVeloClasses	
				if (attsInVeloClassSf(iVc,1,1,iProfile) > 0) then 												 
      				avgAttsInVeloClassSf(iVc,1,iAtt,iProfile) = attsInVeloClassSf(iVc,1,iAtt,iProfile) &
														 / attsInVeloClassSf(iVc,1,1,iProfile) ! 1=nParticles
				end if
			end do
			
		end do ! iAtt
    end do ! iProfile
    
	! Second: get the particle number in units of number of particles L-1 (if in the water 
	! column) and in units of number of particles (if at seafloor)

    do iProfile = 1, nProfiles
    
		!---------------------------------------------------------------------------------
		! Water column recording
		!---------------------------------------------------------------------------------
    
    	do iDh = 1, nImagingDeployDepths
    		
      		do iSc = 1, nSizeClasses
      			if (attsInSizeClass(iSc,iDh,1,iProfile) > 0) then
      				avgAttsInSizeClass(iSc,iDh,1,iProfile) = &
      					(attsInSizeClass(iSc,iDh,1,iProfile)/nAveragingTimeSteps)/cameraSampledVolume(iDh,iProfile) ! # part. L-1   					
      			end if
      		end do
      			
      		do iVc = 1, nVeloClasses
      			if (attsInVeloClass(iVc,iDh,1,iProfile) > 0) then
      				avgAttsInVeloClass(iVc,iDh,1,iProfile) = &
      					(attsInVeloClass(iVc,iDh,1,iProfile)/nAveragingTimeSteps)/cameraSampledVolume(iDh,iProfile) 
      			end if
      		end do
      				
      	end do ! iDh

		!---------------------------------------------------------------------------------
		! Seafloor recording
		!---------------------------------------------------------------------------------
	      	
      	do iSc = 1, nSizeClasses
      		if (attsInSizeClassSf(iSc,1,1,iProfile) > 0) then
      			avgAttsInSizeClassSf(iSc,1,1,iProfile) = attsInSizeClassSf(iSc,1,1,iProfile)/nAveragingTimeSteps       				
      		end if
      	end do
      	
		do iVc = 1, nVeloClasses
			if (attsInVeloClassSf(iVc,1,1,iProfile) > 0) then
      			avgAttsInVeloClassSf(iVc,1,1,iProfile) = attsInVeloClassSf(iVc,1,1,iProfile)/nAveragingTimeSteps
			end if
		end do      	
      	
    end do ! iProfile
    
#ifdef BLOCK_PRINT_INFO_TO_THE_SCREEN
!SPK This really should depend on iProfileMonitor
	
	if (localSeafloorDepth > 1100d0) then
		
		! --- totals ---
		nParticles20m  = SUM(avgAttsInMainType(:,2,1,iProfileMonitor))
		nParticles50m  = SUM(avgAttsInMainType(:,5,1,iProfileMonitor))
		nParticles100m = SUM(avgAttsInMainType(:,10,1,iProfileMonitor))
		nParticles500m = SUM(avgAttsInMainType(:,14,1,iProfileMonitor))
		nParticles1km  = SUM(avgAttsInMainType(:,19,1,iProfileMonitor))
		
		mass20m  = SUM(avgAttsInMainType(:,2,2,iProfileMonitor))
		mass50m  = SUM(avgAttsInMainType(:,5,2,iProfileMonitor))
		mass100m = SUM(avgAttsInMainType(:,10,2,iProfileMonitor))
		mass500m = SUM(avgAttsInMainType(:,14,2,iProfileMonitor))
		mass1km  = SUM(avgAttsInMainType(:,19,2,iProfileMonitor))
		
		fmt2 = '(2X,A,8ES10.2)'
		
		! =====================================================================
		! Fraction of particles per particle type
		! =====================================================================
		
		print *
		print *, 'Fraction of particles per particle type'
		print *, '            Faecal  -  aggs  -  phyto  -  deadphyt -  TEP  -  clay  - zoodead  -  else'
		
		! ---- 20 m ----
		if (nParticles20m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,2,1,iProfileMonitor) / nParticles20m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@20m  :', fracTypes
		
		! ---- 50 m ----
		if (nParticles50m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,5,1,iProfileMonitor) / nParticles50m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@50m  :', fracTypes
		
		! ---- 100 m ----
		if (nParticles100m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,10,1,iProfileMonitor) / nParticles100m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@100m :', fracTypes
		
		! ---- 500 m ----
		if (nParticles500m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,14,1,iProfileMonitor) / nParticles500m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@500m :', fracTypes
		
		! ---- 1 km ----
		if (nParticles1km > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,19,1,iProfileMonitor) / nParticles1km
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@1km  :', fracTypes
		
		! =====================================================================
		! Fraction of mass per particle type
		! =====================================================================
		
		print *
		print *, 'Fraction of mass per particle type'
		print *, '            Faecal  -  aggs  -  phyto  -  deadphyt -  TEP  -  clay  - zoodead  -  else'
		
		! ---- 20 m ----
		if (mass20m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,2,2,iProfileMonitor) / mass20m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@20m  :', fracTypes
		
		! ---- 50 m ----
		if (mass50m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,5,2,iProfileMonitor) / mass50m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@50m  :', fracTypes
		
		! ---- 100 m ----
		if (mass100m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,10,2,iProfileMonitor) / mass100m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@100m :', fracTypes
		
		! ---- 500 m ----
		if (mass500m > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,14,2,iProfileMonitor) / mass500m
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@500m :', fracTypes
		
		! ---- 1 km ----
		if (mass1km > 0d0) then
			fracTypes(:) = avgAttsInMainType(:,19,2,iProfileMonitor) / mass1km
		else
			fracTypes(:) = 0d0
		end if
		print fmt2, '@1km  :', fracTypes
		
		! =====================================================================
		! Mean particle attributes
		! =====================================================================
		
		print *
		print *, 'Average looking particle types: diameter (um)'
		print *, '            Faecal  -  aggs  -  phyto  -  deadphyt -  TEP  -  clay  - zoodead  -  else'
		print fmt2, '@100m :', avgAttsInMainType(:,10,4,iProfileMonitor)
		print fmt2, '@500m :', avgAttsInMainType(:,14,4,iProfileMonitor)
		print fmt2, '@1km  :', avgAttsInMainType(:,19,4,iProfileMonitor)
		
		print *
		print *, 'Average looking particle types: velocity (m d-1)'
		print *, '            Faecal  -  aggs  -  phyto  -  deadphyt -  TEP  -  clay  - zoodead  -  else'
		print fmt2, '@100m :', avgAttsInMainType(:,10,6,iProfileMonitor)
		print fmt2, '@500m :', avgAttsInMainType(:,14,6,iProfileMonitor)
		print fmt2, '@1km  :', avgAttsInMainType(:,19,6,iProfileMonitor)
		
		print *
		print *, 'Average looking particle types: porosity'
		print *, '            Faecal  -  aggs  -  phyto  -  deadphyt -  TEP  -  clay  - zoodead  -  else'
		print fmt2, '@100m :', avgAttsInMainType(:,10,7,iProfileMonitor)
		print fmt2, '@500m :', avgAttsInMainType(:,14,7,iProfileMonitor)
		print fmt2, '@1km  :', avgAttsInMainType(:,19,7,iProfileMonitor)
		
		print *
	end if
	
#endif
    
!#ifndef BLOCK_TMM_SLAMS
!#ifndef SAVE_DISK_SPACE
	newAvgAttsSlice = iLastSliceAvgAtt + 1
	call write_r8_field(nSizeClasses,nImagingDeployDepths*14*nProfiles,newAvgAttsSlice,avgAttsInSizeClass,filenameAvgAttSizeClass) 			
	call write_r8_field(nVeloClasses,nImagingDeployDepths*14*nProfiles,newAvgAttsSlice,avgAttsInVeloClass,filenameAvgAttVeloClass)
	call write_r8_field(nSizeClasses,1*14*nProfiles,newAvgAttsSlice,avgAttsInSizeClassSf,filenameAvgAttSizeClassSf) 			
	call write_r8_field(nVeloClasses,1*14*nProfiles,newAvgAttsSlice,avgAttsInVeloClassSf,filenameAvgAttVeloClassSf)
	call write_r8_field(nMainParticleTypes,nImagingDeployDepths*7*nProfiles,newAvgAttsSlice,avgAttsInMainType,filenameAvgAttMainType)
	call write_r8_field(nMainParticleTypes,1*7*nProfiles,newAvgAttsSlice,avgAttsInMainTypeSf,filenameAvgAttMainTypeSf)		
	iLastSliceAvgAtt = newAvgAttsSlice
!#endif
!#endif	

! 	Reset attribute classifiers
	attsInSizeClass(:,:,:,:)   = 0d0
	attsInVeloClass(:,:,:,:)   = 0d0
	attsInSizeClassSf(:,:,:,:) = 0d0
	attsInVeloClassSf(:,:,:,:) = 0d0
	attsInMainType(:,:,:,:)    = 0d0
	attsInMainTypeSf(:,:,:,:)  = 0d0

end subroutine WriteParticleAverageAttributes

! ========================================================================================

subroutine WriteParticleFlux(localSeafloorDepth,nProfiles)

	!-------------------------------------------------------------------------------------
	! The sediment trap contents have units of mg m-2 collected over a period of time, which
	! is around 30 days. We have to normalise the collection per day. The collection of
	! material happens in the subroutine MarineSnowCatcher, in sink.F90
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nProfiles
	real*8, intent(in) :: localSeafloorDepth
	
	integer :: newSedTrapSlice, iMaterial, i
	character(len=30) :: fmt, fmt2

	do iMaterial = 1, 5 ! 1=orgC, 2=TEP-C, 3=calcite, 4=opal, 5=clay
		particleFlux(iMaterial,:,:) = sedTrap(iMaterial,:,:)/particleCollectionPeriod
		particleFluxSf(iMaterial,1,:) = sedTrapSf(iMaterial,1,:)/particleCollectionPeriod
	end do

#ifndef BLOCK_TMM_SLAMS
	newSedTrapSlice = iLastSliceSedTrap + 1 
 	call write_r8_field(5,nSedTrapDeployDepths*nProfiles,newSedTrapSlice,particleFlux,filenameFlux)
 	call write_r8_field(5,1*nProfiles,newSedTrapSlice,particleFluxSf,filenameFluxSf)
 	iLastSliceSedTrap = newSedTrapSlice
#endif

#ifdef BLOCK_PRINT_INFO_TO_THE_SCREEN
!SPK This really should depend on iProfileMonitor
    
	if (localSeafloorDepth > 1000d0) then	 	
		fmt = '(2X,A,5ES10.1,ES10.1)'	
		fmt2 = '(2X,2ES10.2)'
		print *			
		print *, ' Flux (mg m-2 d-1)'		 
		print *, '          orgC   -  tepC   -   calc   -    opal   -  CaCO3/POC (molar)'
		print *, ' exp:     ~70                  ~100         ~15         0.05-2'
		print fmt, '@100m  :', (particleFlux(i,10,iProfileMonitor), i=1,4), &
			(particleFlux(3,10,iProfileMonitor)/MOLAR_MASS_CACO3)/ &
				((particleFlux(1,10,iProfileMonitor)+particleFlux(2,10,iProfileMonitor))/MOLAR_MASS_CARBON)
		print fmt, '@500m  :', (particleFlux(i,14,iProfileMonitor), i=1,4), &
			(particleFlux(3,14,iProfileMonitor)/MOLAR_MASS_CACO3)/ &
				((particleFlux(1,14,iProfileMonitor)+particleFlux(2,14,iProfileMonitor))/MOLAR_MASS_CARBON)		
		print fmt, '@1km   :', (particleFlux(i,19,iProfileMonitor), i=1,4), &
			(particleFlux(3,19,iProfileMonitor)/MOLAR_MASS_CACO3)/ &
				((particleFlux(1,19,iProfileMonitor)+particleFlux(2,19,iProfileMonitor))/MOLAR_MASS_CARBON)				
		if (localSeafloorDepth > 2000) then	 
		print fmt, '@2km   :', (particleFlux(i,23,iProfileMonitor), i=1,4), &
			(particleFlux(3,23,iProfileMonitor)/MOLAR_MASS_CACO3)/ &
				((particleFlux(1,23,iProfileMonitor)+particleFlux(2,23,iProfileMonitor))/MOLAR_MASS_CARBON)
		if (localSeafloorDepth > 3000) then	
		print fmt, '@3km   :', (particleFlux(i,27,iProfileMonitor), i=1,4), &
			(particleFlux(3,27,iProfileMonitor)/MOLAR_MASS_CACO3)/ &
				((particleFlux(1,27,iProfileMonitor)+particleFlux(2,27,iProfileMonitor))/MOLAR_MASS_CARBON)			
		end if
		end if
		print fmt, '@sf    :', (particleFluxSf(i,1,iProfileMonitor), i=1,4), &
			(particleFluxSf(3,1,iProfileMonitor)/MOLAR_MASS_CACO3)/ &
				((particleFluxSf(1,1,iProfileMonitor)+particleFluxSf(2,1,iProfileMonitor))/MOLAR_MASS_CARBON)
		print *	
		if (localSeafloorDepth > 2000d0) then	
		print *, ' Teff(100-1000m) - Teff(1000-2000m)'
		print fmt2, (particleFlux(1,19,iProfileMonitor)+particleFlux(2,19,iProfileMonitor))/&
				 	(particleFlux(1,10,iProfileMonitor)+particleFlux(2,10,iProfileMonitor)), &
				 	(particleFlux(1,23,iProfileMonitor)+particleFlux(2,23,iProfileMonitor))/&
				 	(particleFlux(1,19,iProfileMonitor)+particleFlux(2,19,iProfileMonitor))
		print *	
		endif		 	
				 	
	end if
#endif

	! Empty the contents of the sediment traps
	sedTrap(:,:,:)   = 0d0
	sedTrapSf(:,:,:) = 0d0

end subroutine WriteParticleFlux

! ========================================================================================

subroutine WriteEulerianVariables(localSeafloorDepth,nProfiles)

	integer, intent(in) :: nProfiles
	real*8, intent(in) :: localSeafloorDepth
	
	real*8 :: waterColumnVolume, waterColumnDepth, molsCtoMgC, &
		nppOrg, nppTepC, nppCaCO3, nppOpal, claySeeded, zIngestInt, zRespInt, mRespInt, &	
		zIngest50m, zEgest50m, zResp50m, zExcret50m, zSolub50m, zDead50m, mResp50m, mSolub50m, &
		zIngest100m, zEgest100m, zResp100m, zExcret100m, zSolub100m, zDead100m, mResp100m, mSolub100m, &
		zIngest500m, zEgest500m, zResp500m, zExcret500m, zSolub500m, zDead500m, mResp500m, mSolub500m, &
		zIngest1km, zEgest1km, zResp1km, zExcret1km, zSolub1km, zDead1km, mResp1km, mSolub1km, &
		zIngest2km, zEgest2km, zResp2km, zExcret2km, zSolub2km, zDead2km, mResp2km, mSolub2km, &
		zooIngestTot, nppTot, picToPoc, tepToPoc, peeff, ratioZooToNpp, totalResp, fracZ
	real*8, dimension(nDepthLayers) :: depthLayerThickness, volume
	integer :: iProfile, iDepthLayer, newEulerianSlice, iSmsTerm	
	character(len=30) :: fmt1, fmt2, fmt3, fmt4, fmt5, fmt6, fmt7, fmt8, fmt9, fmt10, fmt11
	
	avgSMSflux(:,:,:) = 0d0
	avgSMSfluxIntegrated(:,:) = 0d0
	
	do iProfile = 1, nProfiles
	
		!  Normalise per day and volume: mol --> mol m-3 d-1
		do iDepthLayer = 1, nDepthLayers	
			depthLayerThickness(iDepthLayer) = zbot(iDepthLayer) - ztop(iDepthLayer) ! m			
			volume(iDepthLayer) = depthLayerThickness(iDepthLayer)*gridCellArea(iProfile) ! m3				
			avgSMSflux(:,iDepthLayer,iProfile) = accumSMS(:,iDepthLayer,iProfile) &
				/ (volume(iDepthLayer)*particleCollectionPeriod) ! mol m-3 d-1				
		end do
		waterColumnVolume = SUM(volume(:))
		waterColumnDepth = SUM(depthLayerThickness(:))
		
		! Accumulate water-column timestep production: mol --> mol m-2 d-1 
		do iSmsTerm = 1, maxNumSmsTerms
			avgSMSfluxIntegrated(iSmsTerm,iProfile) = &
				(SUM(accumSMS(iSmsTerm,:,iProfile))/(waterColumnVolume*particleCollectionPeriod))*waterColumnDepth ! mol m-2 d-1
		end do
		
	end do
		
!#ifdef BLOCK_PRINT_INFO_TO_THE_SCREEN_EXTENDED

!SPK This really should depend on iProfileMonitor
    molsCtoMgC = MOLAR_MASS_CARBON*1d3
    
	fmt1 = '(2X,A,1ES10.1)'
	fmt2 = '(2X,A,2ES10.1)' 
	fmt3 = '(2X,A,3ES10.1)'
	fmt4 = '(2X,A,4ES10.1)'
	fmt5 = '(2X,A,5ES10.1)'	
	fmt6 = '(2X,A,7ES10.1)'	
	fmt7 = '(2X,2ES10.2)'
	fmt8 = '(2X,3ES10.2)'
	fmt9 = '(2X,4ES10.2)'
	fmt10 = '(2X,5ES10.2)'
	fmt11 = '(2X,6ES10.2)'

	nppOrg = avgSMSfluxIntegrated(iPrimProdOrgC,iProfileMonitor) ! mol m-2 d-1
	nppTepC = avgSMSfluxIntegrated(iProdTepPhyto,iProfileMonitor) ! mol m-2 d-1
	nppCaCO3 = avgSMSfluxIntegrated(iPrimProdCaCO3,iProfileMonitor) ! mol m-2 d-1
	nppOpal = avgSMSfluxIntegrated(iPrimProdOpal,iProfileMonitor) ! mol m-2 d-1
	claySeeded = avgSMSfluxIntegrated(iDepoClay,iProfileMonitor) ! mol m-2 d-1
	
	zooIngestTot = SUM(avgSMSflux(iZooIngestOrgC,:,iProfileMonitor)) + &
               	   SUM(avgSMSflux(iZooIngestTepC,:,iProfileMonitor))
               	   
    nppTot = SUM(avgSMSflux(iPrimProdOrgC,:,iProfileMonitor))
    
    if (ABS(nppOrg) > detection_limit_poc) then
		picToPoc = nppCaCO3 / nppOrg
		tepToPoc = nppTepC  / nppOrg
		peeff = (particleFlux(1,10,iProfileMonitor)+particleFlux(2,10,iProfileMonitor))/(nppOrg*molsCtoMgC)
	else
		picToPoc = 0d0      ! or NaN, or -9999
		tepToPoc = 0d0
		peeff = 0d0
	end if
	
	if (ABS(nppTot) > detection_limit_poc) then
    	ratioZooToNpp = zooIngestTot / nppTot
	else
    	ratioZooToNpp = 0d0    ! or NaN / -9999d0, depending on your convention
	end if
		
	write(*,*)
	write(*,'(a)') 'Seeded material for an avg. day of that month, zeu integrated (mg m-2 d-1)'
	write(*,'(a)') '  orgC      CaCO3      opal       clay    PIC/POC   tepC/orgC'
	write(*,'(2x,*(es10.2))') &
    	nppOrg*molsCtoMgC, &
    	nppCaCO3*MOLAR_MASS_CACO3*1d3, &
    	nppOpal*MOLAR_MASS_OPAL*1d3, &
    	claySeeded*MOLAR_MASS_CLAY*1d3, &
    	picToPoc, &
    	tepToPoc
	write(*,*)
	write(*,'(a45,es12.2)') 'Zoo. ingestion (all water col) / NPP', ratioZooToNpp
	if (localSeafloorDepth > 100d0) write(*,'(a45,es12.2)') 'P_eff (POC/NPP at 100 m)', peeff
	
	if (localSeafloorDepth > 1000d0) then	
	
	zIngestInt = 0d0
	zRespInt = 0d0
	mRespInt = 0d0
	do iDepthLayer = 1, 100
		zIngestInt = zIngestInt + (avgSMSflux(iZooIngestOrgC,iDepthLayer,iProfileMonitor) &
			+ avgSMSflux(iZooIngestTepC,iDepthLayer,iProfileMonitor))*depthLayerThickness(iDepthLayer) ! mol m-2 d-1
		zRespInt = zRespInt + (avgSMSflux(iZooRespOrgC,iDepthLayer,iProfileMonitor) &
			+ avgSMSflux(iZooRespTepC,iDepthLayer,iProfileMonitor))*depthLayerThickness(iDepthLayer) ! mol m-2 d-1
		mRespInt = mRespInt + (avgSMSflux(iMicrobRespOrgC,iDepthLayer,iProfileMonitor) &
			+ avgSMSflux(iMicrobRespTepC,iDepthLayer,iProfileMonitor))*depthLayerThickness(iDepthLayer) ! mol m-2 d-1
	end do
	
	totalResp = zRespInt+mRespInt
	if (abs(totalResp) > detection_limit_poc) then        ! treat tiny/zero as zero
   		fracZ = zRespInt / totalResp
	else
   		fracZ = 0d0                          ! or a fill value like -9999d0
	end if
	fracZ = max(0d0, MIN(1d0, fracZ))
	
	print *	 		
	print *, ' Heterotrophic flows (mg m-2 d-1), integrated over 1000 m'	
	print *, '     zooIng  -  zooResp  -  bacResp  -  zooResp/allResp'
	print fmt9, zIngestInt*molsCtoMgC, zRespInt*molsCtoMgC, mRespInt*molsCtoMgC, fracZ
	
	end if

	if (localSeafloorDepth > 2000d0) then	 
	
	zIngest50m = (avgSMSflux(iZooIngestOrgC,5,iProfileMonitor)+avgSMSflux(iZooIngestTepC,5,iProfileMonitor))*molsCtoMgC
	zEgest50m = (avgSMSflux(iZooEgestOrgC,5,iProfileMonitor)+avgSMSflux(iZooEgestTepC,5,iProfileMonitor))*molsCtoMgC
	zResp50m = (avgSMSflux(iZooRespOrgC,5,iProfileMonitor)+avgSMSflux(iZooRespTepC,5,iProfileMonitor))*molsCtoMgC
	zExcret50m = (avgSMSflux(iZooExcretOrgC,5,iProfileMonitor)+avgSMSflux(iZooExcretTepC,5,iProfileMonitor))*molsCtoMgC
	zSolub50m = (avgSMSflux(iZooSolubOrgC,5,iProfileMonitor)+avgSMSflux(iZooSolubTepC,5,iProfileMonitor))*molsCtoMgC
	zDead50m = avgSMSflux(iZooDeathOrgC,5,iProfileMonitor)*molsCtoMgC
	mResp50m = (avgSMSflux(iMicrobRespOrgC,5,iProfileMonitor)+avgSMSflux(iMicrobRespTepC,5,iProfileMonitor))*molsCtoMgC
	mSolub50m = (avgSMSflux(iMicrobSolubOrgC,5,iProfileMonitor)+avgSMSflux(iMicrobSolubTepC,5,iProfileMonitor))*molsCtoMgC
				
	zIngest100m = (avgSMSflux(iZooIngestOrgC,10,iProfileMonitor)+avgSMSflux(iZooIngestTepC,10,iProfileMonitor))*molsCtoMgC
	zEgest100m = (avgSMSflux(iZooEgestOrgC,10,iProfileMonitor)+avgSMSflux(iZooEgestTepC,10,iProfileMonitor))*molsCtoMgC
	zResp100m = (avgSMSflux(iZooRespOrgC,10,iProfileMonitor)+avgSMSflux(iZooRespTepC,10,iProfileMonitor))*molsCtoMgC
	zExcret100m = (avgSMSflux(iZooExcretOrgC,10,iProfileMonitor)+avgSMSflux(iZooExcretTepC,10,iProfileMonitor))*molsCtoMgC
	zSolub100m = (avgSMSflux(iZooSolubOrgC,10,iProfileMonitor)+avgSMSflux(iZooSolubTepC,10,iProfileMonitor))*molsCtoMgC	
	zDead100m = avgSMSflux(iZooDeathOrgC,10,iProfileMonitor)*molsCtoMgC
	mResp100m = (avgSMSflux(iMicrobRespOrgC,10,iProfileMonitor)+avgSMSflux(iMicrobRespTepC,10,iProfileMonitor))*molsCtoMgC
	mSolub100m = (avgSMSflux(iMicrobSolubOrgC,10,iProfileMonitor)+avgSMSflux(iMicrobSolubTepC,10,iProfileMonitor))*molsCtoMgC

	zIngest500m = (avgSMSflux(iZooIngestOrgC,50,iProfileMonitor)+avgSMSflux(iZooIngestTepC,50,iProfileMonitor))*molsCtoMgC
	zEgest500m = (avgSMSflux(iZooEgestOrgC,50,iProfileMonitor)+avgSMSflux(iZooEgestTepC,50,iProfileMonitor))*molsCtoMgC
	zResp500m = (avgSMSflux(iZooRespOrgC,50,iProfileMonitor)+avgSMSflux(iZooRespTepC,50,iProfileMonitor))*molsCtoMgC
	zExcret500m = (avgSMSflux(iZooExcretOrgC,50,iProfileMonitor)+avgSMSflux(iZooExcretTepC,50,iProfileMonitor))*molsCtoMgC
	zSolub500m = (avgSMSflux(iZooSolubOrgC,50,iProfileMonitor)+avgSMSflux(iZooSolubTepC,50,iProfileMonitor))*molsCtoMgC
	zDead500m = avgSMSflux(iZooDeathOrgC,50,iProfileMonitor)*molsCtoMgC
	mResp500m = (avgSMSflux(iMicrobRespOrgC,50,iProfileMonitor)+avgSMSflux(iMicrobRespTepC,50,iProfileMonitor))*molsCtoMgC
	mSolub500m = (avgSMSflux(iMicrobSolubOrgC,50,iProfileMonitor)+avgSMSflux(iMicrobSolubTepC,50,iProfileMonitor))*molsCtoMgC
	
	zIngest1km = (avgSMSflux(iZooIngestOrgC,100,iProfileMonitor)+avgSMSflux(iZooIngestTepC,100,iProfileMonitor))*molsCtoMgC
	zEgest1km = (avgSMSflux(iZooEgestOrgC,100,iProfileMonitor)+avgSMSflux(iZooEgestTepC,100,iProfileMonitor))*molsCtoMgC
	zResp1km = (avgSMSflux(iZooRespOrgC,100,iProfileMonitor)+avgSMSflux(iZooRespTepC,100,iProfileMonitor))*molsCtoMgC
	zExcret1km = (avgSMSflux(iZooExcretOrgC,100,iProfileMonitor)+avgSMSflux(iZooExcretTepC,100,iProfileMonitor))*molsCtoMgC
	zSolub1km = (avgSMSflux(iZooSolubOrgC,100,iProfileMonitor)+avgSMSflux(iZooSolubTepC,100,iProfileMonitor))*molsCtoMgC		
	zDead1km = avgSMSflux(iZooDeathOrgC,100,iProfileMonitor)*molsCtoMgC
	mResp1km = (avgSMSflux(iMicrobRespOrgC,100,iProfileMonitor)+avgSMSflux(iMicrobRespTepC,100,iProfileMonitor))*molsCtoMgC
	mSolub1km = (avgSMSflux(iMicrobSolubOrgC,100,iProfileMonitor)+avgSMSflux(iMicrobSolubTepC,100,iProfileMonitor))*molsCtoMgC
	
	zIngest2km = (avgSMSflux(iZooIngestOrgC,200,iProfileMonitor)+avgSMSflux(iZooIngestTepC,200,iProfileMonitor))*molsCtoMgC
	zEgest2km = (avgSMSflux(iZooEgestOrgC,200,iProfileMonitor)+avgSMSflux(iZooEgestTepC,200,iProfileMonitor))*molsCtoMgC
	zResp2km = (avgSMSflux(iZooRespOrgC,200,iProfileMonitor)+avgSMSflux(iZooRespTepC,200,iProfileMonitor))*molsCtoMgC
	zExcret2km = (avgSMSflux(iZooExcretOrgC,200,iProfileMonitor)+avgSMSflux(iZooExcretTepC,200,iProfileMonitor))*molsCtoMgC
	zSolub2km = (avgSMSflux(iZooSolubOrgC,200,iProfileMonitor)+avgSMSflux(iZooSolubTepC,200,iProfileMonitor))*molsCtoMgC	
	zDead2km = avgSMSflux(iZooDeathOrgC,200,iProfileMonitor)*molsCtoMgC
	mResp2km = (avgSMSflux(iMicrobRespOrgC,200,iProfileMonitor)+avgSMSflux(iMicrobRespTepC,200,iProfileMonitor))*molsCtoMgC	
	mSolub2km = (avgSMSflux(iMicrobSolubOrgC,200,iProfileMonitor)+avgSMSflux(iMicrobSolubTepC,200,iProfileMonitor))*molsCtoMgC
	
	print *		
	print *, ' Heterotrophic rates (mg m-3 d-1)'
	print *, '            zIng  -  zEgest - zExcret - zMess -   zResp -  mSolub -  mResp'  									  
	print fmt6, '@50m :', zIngest50m, zEgest50m, zExcret50m, zSolub50m, zResp50m, mSolub50m, mResp50m				   		 		 
	print fmt6, '@100m:', zIngest100m, zEgest100m, zExcret100m, zSolub100m, zResp100m, mSolub100m, mResp100m				  	 
	print fmt6, '@500m:', zIngest500m, zEgest500m, zExcret500m, zSolub500m, zResp500m, mSolub500m, mResp500m
	print fmt6, '@1km :', zIngest1km, zEgest1km, zExcret1km, zSolub1km, zResp1km, mSolub1km, mResp1km							  
	print fmt6, '@2km :', zIngest2km, zEgest2km, zExcret2km, zSolub2km, zResp2km, mSolub2km, mResp2km	
	print *		
	print *, ' Fractions '
	print *, '        zS/(zS+R) - bR/(bR+S) - allR/(allS+R)'
	print fmt3, '@50m  :', zResp50m/(zExcret50m+zSolub50m+zResp50m), mResp50m/(mSolub50m+mResp50m), &
						   (zResp50m+mResp50m)/(zExcret50m+zSolub50m+mSolub50m+zResp50m+mResp50m)
	print fmt3, '@100m :', zResp100m/(zExcret100m+zSolub100m+zResp100m), mResp100m/(mSolub100m+mResp100m), &
						   (zResp100m+mResp100m)/(zExcret100m+zSolub100m+mSolub100m+zResp100m+mResp100m)			
	print fmt3, '@500m :', zResp500m/(zExcret500m+zSolub500m+zResp500m), mResp500m/(mSolub500m+mResp500m), &
						   (zResp500m+mResp500m)/(zExcret500m+zSolub500m+mSolub500m+zResp500m+mResp500m)	
	print fmt3, '@1km  :', zResp1km/(zExcret1km+zSolub1km+zResp1km), mResp1km/(mSolub1km+mResp1km), &
						   (zResp1km+mResp1km)/(zExcret1km+zSolub1km+mSolub1km+zResp1km+mResp1km)								   		
	print fmt3, '@2km  :', zResp2km/(zExcret2km+zSolub2km+zResp2km), mResp2km/(mSolub2km+mResp2km), &
						   (zResp2km+mResp2km)/(zExcret2km+zSolub2km+mSolub2km+zResp2km+mResp2km)							   			
	print *	
	print *, '       zR/(zI+zR) - zEx/(zI+zEx) - zEg/(zI+zEg) - zR/(mR+zR)' 							 				 
	print fmt4, '@50m :', zResp50m/(zIngest50m+zResp50m), zExcret50m/(zIngest50m+zExcret50m), &
						  zEgest50m/(zIngest50m+zEgest50m), zResp50m/(mResp50m+zResp50m)			
	print fmt4, '@100m:', zResp100m/(zIngest100m+zResp100m), zExcret100m/(zIngest100m+zExcret100m), &
					      zEgest100m/(zIngest100m+zEgest100m), zResp100m/(mResp100m+zResp100m)						 		
	print fmt4, '@10m :', zResp500m/(zIngest500m+zResp500m), zExcret500m/(zIngest500m+zExcret500m), &
						  zEgest500m/(zIngest500m+zEgest500m), zResp500m/(mResp500m+zResp500m)	
	print fmt4, '@1km :', zResp1km/(zIngest1km+zResp1km), zExcret1km/(zIngest1km+zExcret1km), &
						  zEgest1km/(zIngest1km+zEgest1km), zResp1km/(mResp1km+zResp1km)						 	
	print fmt4, '@2km :', zResp2km/(zIngest2km+zResp2km), zExcret2km/(zIngest2km+zExcret2km), &
						  zEgest2km/(zIngest2km+zEgest2km), zResp2km/(mResp2km+zResp2km)		 
	print *
	print *, ' Mineral dissolution (mg m-3 d-1)'
	print *, '          zooDisCalc - abioDisCalc - abioDisOpal' 
	print fmt3, '@10m :', avgSMSflux(iZooDissolCaCO3,1,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolCaCO3,1,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolOpal,1,iProfileMonitor)*MOLAR_MASS_OPAL*1d3					  							  					  
	print fmt3, '@100m:', avgSMSflux(iZooDissolCaCO3,10,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolCaCO3,10,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolOpal,10,iProfileMonitor)*MOLAR_MASS_OPAL*1d3			
	print fmt3, '@1km :', avgSMSflux(iZooDissolCaCO3,100,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolCaCO3,100,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolOpal,100,iProfileMonitor)*MOLAR_MASS_OPAL*1d3			
	print fmt3, '@2km :', avgSMSflux(iZooDissolCaCO3,200,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolCaCO3,200,iProfileMonitor)*MOLAR_MASS_CACO3*1d3, &
						  avgSMSflux(iDissolOpal,200,iProfileMonitor)*MOLAR_MASS_OPAL*1d3	

	end if

!#endif

#ifndef BLOCK_TMM_SLAMS
!#ifndef SAVE_DISK_SPACE
	newEulerianSlice = iLastSliceEulerian + 1
	call write_r8_field(maxNumSmsTerms,nDepthLayers*nProfiles,newEulerianSlice,avgSMSflux,filenameSms) ! mol m-3 d-1
	call write_r8_field(maxNumSmsTerms,nProfiles,newEulerianSlice,avgSMSfluxIntegrated,filenameSmsIntegrated) ! mol m-2 d-1
	iLastSliceEulerian = newEulerianSlice 
!#endif
#endif

	! Reset values to zero
	accumSMS(:,:,:) = 0d0

end subroutine WriteEulerianVariables

! ========================================================================================

subroutine WriteAuxiliaryTerms(localSeafloorDepth,nProfiles)

	integer, intent(in) :: nProfiles
	real*8, intent(in) :: localSeafloorDepth
	
	integer :: iProfile, iAuxTerm, iDepthLayer, newAuxTermsSlice
	real*8, dimension(:,:,:), allocatable :: avgAuxTerm
	real*8 :: molsCtoMgC, totalPhytoBiomass
	character(len=30) :: fmt1, fmt2, fmt3
	
	allocate(avgAuxTerm(maxNumAuxTerms,nDepthLayers,nProfiles))
	avgAuxTerm(:,:,:) = 0d0

	do iProfile = 1, nProfiles
		do iDepthLayer = 1, nDepthLayers			
			do iAuxTerm = 1, maxNumAuxTerms		
				if (auxCount(iAuxTerm,iDepthLayer,iProfile) > 0) then
				
					! For the zooplankton encounter stats
					if (iAuxTerm == iNumParticlesEvalEncounter) then 
						avgAuxTerm(iAuxTerm,iDepthLayer,iProfile) = auxCount(iAuxTerm,iDepthLayer,iProfile)
					else if (iAuxTerm == iNumParticlesEncountered) then
						avgAuxTerm(iAuxTerm,iDepthLayer,iProfile) = auxCount(iNumParticlesEncountered,iDepthLayer,iProfile)&
							/auxCount(iNumParticlesEvalEncounter,iDepthLayer,iProfile)
					else if (iAuxTerm > iNumParticlesEncountered) then
						avgAuxTerm(iAuxTerm,iDepthLayer,iProfile) = auxCount(iAuxTerm,iDepthLayer,iProfile)&
							/auxCount(iNumParticlesEncountered,iDepthLayer,iProfile)
				
					! For the rest that is not zooplankton (uses default units of the aux term in question)				
					else		
						avgAuxTerm(iAuxTerm,iDepthLayer,iProfile) = auxTerm(iAuxTerm,iDepthLayer,iProfile) &
							/ auxCount(iAuxTerm,iDepthLayer,iProfile) 
					end if
											
				end if
			end do ! iAuxTerm
		end do ! iDepthLayer
	end do ! iProfile

	if (localSeafloorDepth > 100d0) then			
		print *, 'Avg. spec. resp. rate @100 m mesozoo (d-1) ....', avgAuxTerm(iZooSpecRespRate,10,iProfileMonitor)*SECONDS_PER_DAY
		print *, '............................ microb (d-1) .....', avgAuxTerm(iMicrobSpecRespRate,10,iProfileMonitor)*SECONDS_PER_DAY
		print *
		print *, 'Frac. particles fragmented ....................', avgAuxTerm(iNumParticlesEncountered+1,10,iProfileMonitor)
		print *, 'Frac. particles ingested ......................', avgAuxTerm(iNumParticlesEncountered+2,10,iProfileMonitor)
		print *, 'Frac. particles omitted .......................', avgAuxTerm(iNumParticlesEncountered+3,10,iProfileMonitor)
		print *
	end if

	fmt1 = '(2X,A,4ES10.1)'
	fmt2 = '(2X,A,2ES10.1)'
	fmt3 = '(2X,2ES10.1)'
	
	molsCtoMgC = MOLAR_MASS_CARBON*1d3

	totalPhytoBiomass = avgAuxTerm(iDiatBiomass,1,iProfileMonitor) &
					  + avgAuxTerm(iFlagelBiomass,1,iProfileMonitor) &
					  + avgAuxTerm(iCoccoBiomass,1,iProfileMonitor) &
					  + avgAuxTerm(iPicoBiomass,1,iProfileMonitor)
					  
	print *
 	print *, 'Phyto biomass (mg C m-3): diatom .............', avgAuxTerm(iDiatBiomass,1,iProfileMonitor)*molsCtoMgC
 	print *, '......................... flagel .............', avgAuxTerm(iFlagelBiomass,1,iProfileMonitor)*molsCtoMgC
 	print *, '......................... cocco ..............', avgAuxTerm(iCoccoBiomass,1,iProfileMonitor)*molsCtoMgC					 
 	print *, '......................... pico ...............', avgAuxTerm(iPicoBiomass,1,iProfileMonitor)*molsCtoMgC
	print *
 	print *, 'Phyto relat prob: diatom .....................', avgAuxTerm(iAvgProbDiat,1,iProfileMonitor)
 	print *, '................. flagel .....................', avgAuxTerm(iAvgProbFlagel,1,iProfileMonitor)
 	print *, '................. cocco ......................', avgAuxTerm(iAvgProbCocco,1,iProfileMonitor)					 
 	print *, '................. pico .......................', avgAuxTerm(iAvgProbPico,1,iProfileMonitor)
 	print *
 	print *, 'Euphotic depth (m) ...........................', avgAuxTerm(iEuphoticDepth,1,iProfileMonitor)
	print *
	print *, 'Zoo. migration depths (m): night upper .......', avgAuxTerm(iNightDvmUpperBound,1,iProfileMonitor)
	print *, '.......................... night lower .......', avgAuxTerm(iNightDvmLowerBound,1,iProfileMonitor)
	print *, '.......................... day upper .........', avgAuxTerm(iDayDvmUpperBound,1,iProfileMonitor)
	print *, '.......................... day lower .........', avgAuxTerm(iDayDvmLowerBound,1,iProfileMonitor)
	print *

#ifdef BLOCK_PRINT_INFO_TO_THE_SCREEN_EXTENDED

	if (localSeafloorDepth > 1000d0) then			 
	 	print *
		print *, 'Mean zoo. biomass (mg C m-3) throughout the interval 0-200 m'
		print fmt3, SUM(avgAuxTerm(iZooBiomass,1:20,iProfileMonitor)*molsCtoMgC)/20d0

	end if
	if (localSeafloorDepth > 100d0) then
		print *
		print *, 'Total number of zoo. individuals (ind. m-3)'
		print fmt2, '@10m :', avgAuxTerm(iZooNumber,1,iProfileMonitor)
		print fmt2, '@100m:', avgAuxTerm(iZooNumber,10,iProfileMonitor)		 
	end if
	  	! print *	  
! 		print *, 'Relative phyto biomass: diatom ...............', avgAuxTerm(iDiatBiomass,1,iProfileMonitor)/totalPhytoBiomass
! 		print *, '....................... flagel ...............', avgAuxTerm(iFlagelBiomass,1,iProfileMonitor)/totalPhytoBiomass
! 		print *, '....................... cocco ................', avgAuxTerm(iCoccoBiomass,1,iProfileMonitor)/totalPhytoBiomass					 
! 		print *, '....................... pico .................', avgAuxTerm(iPicoBiomass,1,iProfileMonitor)/totalPhytoBiomass						 	  
! 		print *
! 		print *, 'Avg. mol C per cell: diatom ..................', avgAuxTerm(iFreshDiatCellQuota,1,iProfileMonitor)
! 		print *, '.................... flagel ...................', avgAuxTerm(iFreshFlagelCellQuota,1,iProfileMonitor)
! 		print *, '.................... cocco ...................', avgAuxTerm(iFreshCoccoCellQuota,1,iProfileMonitor)
! 		print *, '.................... pico ....................', avgAuxTerm(iFreshPicoCellQuota,1,iProfileMonitor)
! 		print *
	if (localSeafloorDepth > 1000d0) then
		print *, 'Coagulation success @ 50 m ...................', avgAuxTerm(iCoagulationSuccess,5,iProfileMonitor)
		print *, '................... @ 100 m ..................', avgAuxTerm(iCoagulationSuccess,10,iProfileMonitor)
		print *, '................... @ 200 m ..................', avgAuxTerm(iCoagulationSuccess,20,iProfileMonitor)
		print *, '................... @ 500 m ..................', avgAuxTerm(iCoagulationSuccess,50,iProfileMonitor)
		print *, '................... @ 1 km ...................', avgAuxTerm(iCoagulationSuccess,100,iProfileMonitor)
		print *	
		print *, 'Avg collision kernel @ 100 m (m3 s-1) ........', avgAuxTerm(iCollisionKernel,10,iProfileMonitor)
		print *, 'Brownian motion kernel .......................', avgAuxTerm(iBrownianKernel,10,iProfileMonitor)
		print *, 'Differential settling kernel .................', avgAuxTerm(iSettlingKernel,10,iProfileMonitor)
		print *, 'Shear stress kernel ..........................', avgAuxTerm(iShearKernel,10,iProfileMonitor)
		print *, 'Avg collision kernel @ 1 km ..................', avgAuxTerm(iCollisionKernel,100,iProfileMonitor)
		print *	
	end if
	
#endif

#ifndef BLOCK_TMM_SLAMS
	newAuxTermsSlice = iLastSliceAux + 1
	call write_r8_field(maxNumAuxTerms,nDepthLayers*nProfiles,newAuxTermsSlice,avgAuxTerm,filenameAuxTerms)
	iLastSliceAux = newAuxTermsSlice
#endif

	deallocate(avgAuxTerm)

	! Reset values to zero
	auxTerm(:,:,:) = 0d0
	auxCount(:,:,:) = 0d0

end subroutine WriteAuxiliaryTerms

! ========================================================================================

subroutine WriteAndPrintModelStatistics(particle, nClusters, iLastLocus, iTimeStep, nProfiles)

	integer, intent(in) :: nClusters, iLastLocus, iTimeStep, nProfiles
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	
	integer :: newStatsSlice, iCluster, iProfile
	real*8 :: nParticles, collisionSuccess, grazingSuccess

	write(*,*) '----------------------------------------------------------------------------'
	write(*,*)
	write(*,*) '-------------------- Time step  ', iTimeStep , '--------------------------'		
	write(*,*) '-------------------- No clusters', iLastLocus, '--------------------------'

	nParticles = 0d0
	do iCluster = 1, iLastLocus
		if (particle(iCluster)%phase == 1) then
			nParticles = nParticles + particle(iCluster)%nPxC
		end if		
	end do

#ifndef BLOCK_TMM_SLAMS	
#ifndef SAVE_DISK_SPACE
	newStatsSlice = iLastSliceStats + 1
	call write_i4_field(1,1,newStatsSlice,iLastLocus,filenameStatsNumClusters)
	call write_r8_field(1,1,newStatsSlice,nParticles,filenameStatsNumParticles)
	iLastSliceStats = newStatsSlice
#endif
#endif

    do iProfile = 1, nProfiles
    
	  if (nTimesEnteringCollLoop(iProfile) > 0) then
      	collisionSuccess = 100d0*(nCollisionsType1(iProfile)+nCollisionsType2(iProfile) &
		  + nCollisionsType3(iProfile)+nCollisionsType4(iProfile))/REAL(nTimesEnteringCollLoop(iProfile))
	  else
		collisionSuccess = 0d0 
	  end if
      if (nParticlesEvaluatedForGrazing(iProfile) > 0) then
   		grazingSuccess = 100d0*nZooIngestionEvents(iProfile)/REAL(nParticlesEvaluatedForGrazing(iProfile))
	  else
   	    grazingSuccess = 0d0
	  end if
		
	  write(*,*) 'Profile# = ',iProfile
	  write(*,*) 'No. particles (billion) ......................', nParticles/1d9
	  write(*,*) 'No. collisions type 1 ........................', nCollisionsType1(iProfile)
	  write(*,*) '.............. type 2 ........................', nCollisionsType2(iProfile)
	  write(*,*) '.............. type 3 ........................', nCollisionsType3(iProfile)
	  write(*,*) '.............. type 4 ........................', nCollisionsType4(iProfile)
	  write(*,*) 'Collision success, % .........................', collisionSuccess
	  write(*,*) 'No. mesozoo ingestion events .................', nZooIngestionEvents(iProfile)
	  write(*,*) 'No. mesozoo fragmentation events .............', nZooFragmentationEvents(iProfile)
	  write(*,*) 'Grazing success, % ...........................', grazingSuccess
	  write(*,*) 'No. faecal clusters produced .................', nFaecalPelletClustersProduced(iProfile)	
	  write(*,*) 'No. dead zoo clusters produced ...............', nZooDeadClustersProduced(iProfile)
	  write(*,*) 'No. big agg clusters fragmented ..............', nFragmentedBigClusters(iProfile)
	  write(*,*) 'No. mineral clusters dissolved ...............', nMineralClustersDissolved(iProfile)
	  write(*,*) 'No. tiny agg clusters respired ...............', nRespiredTinyClusters(iProfile)
	  write(*,*) 'No. clusters remineralised by microbes .......', nMicrobRespiredClusters(iProfile)
	  write(*,*) 'No. clusters photolysed ......................', nClustersPhotolysed(iProfile)
	  write(*,*) 'No. clusters that have reached seafloor ......', nClustersAtSeafloor(iProfile)

    end do ! iProfile
  
end subroutine WriteAndPrintModelStatistics

! ========================================================================================

subroutine WriteSnapshotsOfParticleAttributes(nPhase1clusters, phase1clusterIndices, particle, &
	nClusters, nProfiles)
	
	integer, intent(in) :: nClusters, nPhase1clusters, nProfiles
	integer, dimension(nPhase1clusters), intent(in) :: phase1clusterIndices
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle

	integer :: iSlice, iCluster, iPhase1Cluster
	integer, dimension(:,:), allocatable :: tempi
	real*8, dimension(:,:), allocatable :: tempr

	iLastSliceSnapshot = iLastSliceSnapshot + 1
	iSlice = iLastSliceSnapshot
	
	allocate(tempi(nClusters,7))
	allocate(tempr(nClusters,19))
	tempi(:,:) = 0
	tempr(:,:) = 0d0
	
	do iPhase1Cluster = 1, nPhase1clusters
		iCluster = phase1clusterIndices(iPhase1Cluster)
	
		tempi(iPhase1Cluster,1) = particle(iCluster)%id							!1
		tempi(iPhase1Cluster,2) = particle(iCluster)%initType					!2
		tempi(iPhase1Cluster,3) = particle(iCluster)%initPft					!3
		tempi(iPhase1Cluster,4) = particle(iCluster)%phase						!4
		tempi(iPhase1Cluster,5) = particle(iCluster)%living						!5
		tempi(iPhase1Cluster,6) = particle(iCluster)%faecal						!6
		tempi(iPhase1Cluster,7) = particle(iCluster)%tstepCreat					!7
		
		tempr(iPhase1Cluster,1) = particle(iCluster)%molesOrgC					!8
		tempr(iPhase1Cluster,2) = particle(iCluster)%molesTepC					!9
		tempr(iPhase1Cluster,3) = particle(iCluster)%molesMineral(iOpal)		!10
		tempr(iPhase1Cluster,4) = particle(iCluster)%molesMineral(iCalcite)		!11
		tempr(iPhase1Cluster,5) = particle(iCluster)%molesMineral(iClay)		!12
		tempr(iPhase1Cluster,6) = particle(iCluster)%massOrgMatter				!13
		tempr(iPhase1Cluster,7) = particle(iCluster)%massTep					!14
		tempr(iPhase1Cluster,8) = particle(iCluster)%nPxC						!15
		tempr(iPhase1Cluster,9) = particle(iCluster)%nPpxP						!16
		tempr(iPhase1Cluster,10) = particle(iCluster)%nPpxC						!17
		tempr(iPhase1Cluster,11) = particle(iCluster)%depth						!18
		tempr(iPhase1Cluster,12) = particle(iCluster)%mass						!19
		tempr(iPhase1Cluster,13) = particle(iCluster)%solidVolume				!20
		tempr(iPhase1Cluster,14) = particle(iCluster)%density					!21
		tempr(iPhase1Cluster,15) = particle(iCluster)%radius					!22
		tempr(iPhase1Cluster,16) = particle(iCluster)%radiusPp					!23
		tempr(iPhase1Cluster,17) = particle(iCluster)%porosity					!24
		tempr(iPhase1Cluster,18) = particle(iCluster)%stickiness				!25
		tempr(iPhase1Cluster,19) = particle(iCluster)%velocity					!26
	end do

!#ifdef SAVE_DISK_SPACE	
	call write_i4_field(nClusters,7*nProfiles,iSlice,tempi,filenameClustersInteger)
	call write_r8_field(nClusters,19*nProfiles,iSlice,tempr,filenameClustersReal)	
!#endif
		
	deallocate(tempi)
	deallocate(tempr)
	
end subroutine WriteSnapshotsOfParticleAttributes

! ========================================================================================

subroutine WriteInstantaneousSnapshots(localSeafloorDepth, nProfiles)

	integer, intent(in) :: nProfiles
	real*8, intent(in) :: localSeafloorDepth
	
	integer :: iDh, iSc, iVc, iAtt, iTy, iSlice, iProfile, i
    character(len=30) :: fmt1, fmt2
    
    avgAttsInSizeClassInst(:,:,:,:)   = 0d0
	avgAttsInVeloClassInst(:,:,:,:)   = 0d0
	avgAttsInSizeClassSfInst(:,:,:,:) = 0d0
	avgAttsInVeloClassSfInst(:,:,:,:) = 0d0
	
	avgAttsInMainTypeInst(:,:,:,:)    = 0d0
	avgAttsInMainTypeSfInst(:,:,:,:)  = 0d0

! 	Slices 2-14 of array avgAtts: get the average attribute for that particle size class 
!	or velo class by dividing the accum attribute by the number of particles

    do iProfile = 1, nProfiles

	!	Water column
		do iDh = 1, nImagingDeployDepths
		
		!	Particle type classifier
			do iTy = 1, nMainParticleTypes
				avgAttsInMainTypeInst(iTy,iDh,1,iProfile) = attsInMainTypeInst(iTy,iDh,1,iProfile) ! number of particles
				avgAttsInMainTypeInst(iTy,iDh,2,iProfile) = attsInMainTypeInst(iTy,iDh,2,iProfile) ! total mass
				avgAttsInMainTypeInst(iTy,iDh,3,iProfile) = attsInMainTypeInst(iTy,iDh,3,iProfile) ! total POC (orgC+TEPC)
				do iAtt = 4, 7 ! avg. diameter of a particle, avg. density of a particle, avg. velocity of a particle, avg. porosity of a particle
					if (attsInMainTypeInst(iTy,iDh,1,iProfile) > 0) then
						avgAttsInMainTypeInst(iTy,iDh,iAtt,iProfile) = attsInMainTypeInst(iTy,iDh,iAtt,iProfile) &
							/attsInMainTypeInst(iTy,iDh,1,iProfile) ! 1=nParticles
					end if				
				end do
			end do		

			do iAtt = 2, 14

		 	! 	Size classifier in the water column	
				do iSc = 1, nSizeClasses
				  	if (attsInSizeClassInst(iSc,iDh,1,iProfile) > 0) then
						avgAttsInSizeClassInst(iSc,iDh,iAtt,iProfile) = attsInSizeClassInst(iSc,iDh,iAtt,iProfile) &
														 / attsInSizeClassInst(iSc,iDh,1,iProfile) ! 1=nParticles
				  	end if				
			  	end do
								
		  	! 	Velocity classifier	in the water column	
			  	do iVc = 1, nVeloClasses	
				  	if (attsInVeloClassInst(iVc,iDh,1,iProfile) > 0) then 	
					  	avgAttsInVeloClassInst(iVc,iDh,iAtt,iProfile) = attsInVeloClassInst(iVc,iDh,iAtt,iProfile) &
														 / attsInVeloClassInst(iVc,iDh,1,iProfile) ! 1=nParticles
				  	end if	
			  	end do

		  	end do ! iAtt	
				
	  	end do ! iDh

	!	Seafloor
	
	!	Particle type classifier	
		do iTy = 1, nMainParticleTypes
			do iAtt = 4, 7 ! avg. diameter of a particle, avg. density of a particle, avg. velocity of a particle, avg. porosity of a particle
				if (attsInMainTypeSfInst(iTy,1,1,iProfile) > 0) then
					avgAttsInMainTypeSfInst(iTy,1,iAtt,iProfile) = attsInMainTypeSfInst(iTy,1,iAtt,iProfile) &
						/attsInMainTypeSfInst(iTy,1,1,iProfile) ! 1=nParticles
				end if				
			end do
			avgAttsInMainTypeSfInst(iTy,1,1,iProfile) = attsInMainTypeSfInst(iTy,1,1,iProfile) ! number of particles
			avgAttsInMainTypeSfInst(iTy,1,2,iProfile) = attsInMainTypeSfInst(iTy,1,2,iProfile) ! total mass
			avgAttsInMainTypeSfInst(iTy,1,3,iProfile) = attsInMainTypeSfInst(iTy,1,3,iProfile) ! total organic carbon
		end do
			
		do iAtt = 2, 14
		
		!	Size classifier at seafloor		
			do iSc = 1, nSizeClasses
				if (attsInSizeClassSfInst(iSc,1,1,iProfile) > 0) then
	  				avgAttsInSizeClassSfInst(iSc,1,iAtt,iProfile) = attsInSizeClassSfInst(iSc,1,iAtt,iProfile) &
														 / attsInSizeClassSfInst(iSc,1,1,iProfile) ! 1=nParticles
				end if
			end do
				
		!	Velocity classifier at seafloor			
			do iVc = 1, nVeloClasses	
				if (attsInVeloClassSfInst(iVc,1,1,iProfile) > 0) then 												 
      				avgAttsInVeloClassSfInst(iVc,1,iAtt,iProfile) = attsInVeloClassSfInst(iVc,1,iAtt,iProfile) &
														 / attsInVeloClassSfInst(iVc,1,1,iProfile) ! 1=nParticles
				end if
			end do
			
		end do ! iAtt
			
    end do ! iProfile
    
! 	Slice 1 of array avgAtts: get the particle number in units of number of particles L-1
!	(if in the water column) and in units of number of particles (if at seafloor)

    do iProfile = 1, nProfiles
    
    !	For the water column samples
    
    	do iDh = 1, nImagingDeployDepths
    		
      		do iSc = 1, nSizeClasses
      			if (attsInSizeClassInst(iSc,iDh,1,iProfile) > 0) then
      				if (iDh == 10) then
      					print *, 'iSc', iSc, 'pNum', attsInSizeClassInst(iSc,iDh,1,iProfile), 'vol in L', cameraSampledVolume(iDh,iProfile)
      				end if
      				avgAttsInSizeClassInst(iSc,iDh,1,iProfile) = &
      					attsInSizeClassInst(iSc,iDh,1,iProfile)/cameraSampledVolume(iDh,iProfile)     					
      			end if
      		end do
      			
      		do iVc = 1, nVeloClasses
      			if (attsInVeloClassInst(iVc,iDh,1,iProfile) > 0) then
      				avgAttsInVeloClassInst(iVc,iDh,1,iProfile) = &
      					attsInVeloClassInst(iVc,iDh,1,iProfile)/cameraSampledVolume(iDh,iProfile) 
      			end if
      		end do
      				
      	end do ! iDh

	!	For the seafloor samples
	      	
      	do iSc = 1, nSizeClasses
      		if (attsInSizeClassSfInst(iSc,1,1,iProfile) > 0) then
      			avgAttsInSizeClassSfInst(iSc,1,1,iProfile) = attsInSizeClassSfInst(iSc,1,1,iProfile)       				
      		end if
      	end do
      	
		do iVc = 1, nVeloClasses
			if (attsInVeloClassSfInst(iVc,1,1,iProfile) > 0) then
      			avgAttsInVeloClassSfInst(iVc,1,1,iProfile) = attsInVeloClassSfInst(iVc,1,1,iProfile)
			end if
		end do      	
      	
    end do ! iProfile
    
!#ifdef BLOCK_PRINT_INFO_TO_THE_SCREEN
!SPK This really should depend on iProfileMonitor
	
	if (localSeafloorDepth > 1100d0) then		
		fmt1 = '(2X,11ES10.2)'
		print *		
		print *, ' INSTANTANEOUS particle num (#/L) in these size (diameter, um) categories @ 100 m:'	
		print *, '      8.0   -   16.0    -   32.0   -  64.0  -  128'
		print *, ' expected #/L:                      ~100'
		print fmt1, (avgAttsInSizeClassInst(i,10,1,iProfileMonitor), i=4,8)
		print *	
		print *, '      256   -   512  -   1024  -   2048  -   4096'
		print fmt1, (avgAttsInSizeClassInst(i,10,1,iProfileMonitor), i=9,13)
		print *
	end if
	
!#endif
    
!#ifndef BLOCK_TMM_SLAMS
!#ifndef SAVE_DISK_SPACE

! 	Write to files
	iLastSliceParticleNumSnapshot = iLastSliceParticleNumSnapshot + 1
	iSlice = iLastSliceParticleNumSnapshot
	
	call write_r8_field(nSizeClasses,nImagingDeployDepths*14*nProfiles,iSlice,avgAttsInSizeClassInst,filenameInstAvgAttSizeClass) 			
	call write_r8_field(nVeloClasses,nImagingDeployDepths*14*nProfiles,iSlice,avgAttsInVeloClassInst,filenameInstAvgAttVeloClass)

	call write_r8_field(nSizeClasses,1*14*nProfiles,iSlice,avgAttsInSizeClassSfInst,filenameInstAvgAttSizeClassSf) 			
	call write_r8_field(nVeloClasses,1*14*nProfiles,iSlice,avgAttsInVeloClassSfInst,filenameInstAvgAttVeloClassSf)
	
	call write_r8_field(nMainParticleTypes,nImagingDeployDepths*7*nProfiles,iSlice,avgAttsInMainTypeInst,filenameInstAvgAttMainType)
	call write_r8_field(nMainParticleTypes,1*7*nProfiles,iSlice,avgAttsInMainTypeSfInst,filenameInstAvgAttMainTypeSf)		

!#endif
!#endif

! 	Reset attribute classifiers
	attsInSizeClassInst(:,:,:,:)   = 0d0
	attsInVeloClassInst(:,:,:,:)   = 0d0
	attsInSizeClassSfInst(:,:,:,:) = 0d0
	attsInVeloClassSfInst(:,:,:,:) = 0d0
	attsInMainTypeInst(:,:,:,:)    = 0d0
	attsInMainTypeSfInst(:,:,:,:)  = 0d0	

end subroutine WriteInstantaneousSnapshots

! ========================================================================================

subroutine GetParticleAverageAttributesByVolumeClass(volumeClass, particle, nClusters, &
	iCluster, iTimeStep)

	!-------------------------------------------------------------------------------------
	! Subroutine only to be used with Jackson's model
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster, iTimeStep
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	real*8, dimension(nVolumeClasses,11), intent(inout) :: volumeClass
		
	integer :: iVolc, iDh
	real*8 :: np
	
	iVolc = 0
	iDh   = 0
	np = dble(particle(iCluster)%nPxC)
	
	iVolc = FindParticleVolumeClass(particle(iCluster)%solidVolume)
	if (iVolc == 0) then ! particle too large
		lossTerms(2,nVolumeClasses) = lossTerms(2,nVolumeClasses) + np
	end if
	if (particle(iCluster)%depth < 50d0) then ! same layer thickness as in Jackson's model
		iDh = 1		
	end if
	
	if (iDh > 0 .and. iVolc > 0) then
		volumeClass(iVolc,1) = volumeClass(iVolc,1) + np
		volumeClass(iVolc,2) = volumeClass(iVolc,2) + particle(iCluster)%density*np
		volumeClass(iVolc,3) = volumeClass(iVolc,3) + particle(iCluster)%velocity*np						 
		volumeClass(iVolc,4) = volumeClass(iVolc,4) + 2d0*particle(iCluster)%radius*np					 			
		volumeClass(iVolc,5) = volumeClass(iVolc,5) + particle(iCluster)%stickiness*np	
		volumeClass(iVolc,6) = volumeClass(iVolc,6) + particle(iCluster)%porosity*np
		volumeClass(iVolc,7) = volumeClass(iVolc,7) + particle(iCluster)%solidVolume*np
		volumeClass(iVolc,8) = volumeClass(iVolc,8) + particle(iCluster)%fracDim*np
		volumeClass(iVolc,9) = volumeClass(iVolc,9) + 2d0*particle(iCluster)%radiusPp*np
		volumeClass(iVolc,10) = volumeClass(iVolc,10) + particle(iCluster)%depth*np				
		volumeClass(iVolc,11) = volumeClass(iVolc,11) + particle(iCluster)%molesOrgC*np	
	end if

end subroutine GetParticleAverageAttributesByVolumeClass

! ========================================================================================

subroutine WriteTestSmoluchowski(nProfiles, iTimeStep)

	integer, intent(in) :: nProfiles, iTimeStep
	
	integer :: newAvgVolSlice, iVolc, iAtt, iProfile
	real*8 :: nAccumParticles, vol
    character(len=30) :: fmt1, fmt
    
    fmt = '(7F9.1)'
	
	avgAttsInVolumeClass(:,:,:) = 0d0

	! ---- Attributes 2–11: averages ----
    do iProfile = 1, nProfiles		
		do iAtt = 2, 11	
			do iVolc = 1, nVolumeClasses
				nAccumParticles = attsInVolumeClass(iVolc,1,iProfile)
				if (nAccumParticles <= 0) cycle
				avgAttsInVolumeClass(iVolc,iAtt,iProfile) = attsInVolumeClass(iVolc,iAtt,iProfile)/nAccumParticles			
			end do
		end do	
    end do
    
	! ---- Attribute 1: number concentration ----
    do iProfile = 1, nProfiles
    	vol = cameraSampledVolume(1,iProfile)
		do iVolc = 1, nVolumeClasses
			nAccumParticles = attsInVolumeClass(iVolc,1,iProfile)
			if (nAccumParticles <= 0) cycle
			avgAttsInVolumeClass(iVolc,1,iProfile) = nAccumParticles/vol ! num. L-1     					
		end do
    end do
    
	! Write to files
	newAvgVolSlice = iLastSliceAvgVolAtt + 1
	call write_r8_field(nVolumeClasses,11*nProfiles,newAvgVolSlice,avgAttsInVolumeClass,filenameAvgAttVolumeClass)
	iLastSliceAvgVolAtt = newAvgVolSlice
	
	! Reset to 0
	attsInVolumeClass(:,:,:) = 0d0

end subroutine WriteTestSmoluchowski

! ========================================================================================

subroutine WriteModelClosureInformation(iTimeStep, iLastCluster)

	integer, intent(in) :: iTimeStep, iLastCluster
	integer, dimension(9) :: control

	write(*,*)
	write(*,*) 'Writing model closure information...'
	write(*,*)
	
	control = (/iTimeStep-1,                 &			
		    	iLastCluster, 				 &				
				iLastSliceSedTrap,           &
				nSedTrapDeployDepths,        &
				nImagingDeployDepths,        &
				iLastSliceAvgAtt,            &				
				iLastSliceStats,             &
				iLastSliceAux,               &
				iLastSliceSnapshot           /)

	call write_i4_field(1,9,1,control,filenameControl)
	
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	call write_r8_field(2,nVolumeClasses,1,lossTerms,filenameLossTerms) 
#endif

end subroutine WriteModelClosureInformation

! ========================================================================================

end module modeloutput