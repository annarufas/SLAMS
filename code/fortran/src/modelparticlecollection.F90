#include "blockdefinitions.h"

module modelparticlecollection

! ----------------------------------------------------------------------------------------
! This module handles the arrays that collect particle flux data and classify the
! particles collected according to specific size and velocity categories.
! ----------------------------------------------------------------------------------------

use modelparameters, only: maxNumSedTrapDeployDepths, maxNumImagingDeployDepths, &
	nDepthLayers, possibleSedTrapDeployDepths, possibleImagingDeployDepths, nSizeClasses, &
	nVeloClasses, nVolumeClasses, nMainParticleTypes
use modelgrid, only: gridCellArea, ztop, zbot

implicit none
public

integer, public :: nSedTrapDeployDepths, nImagingDeployDepths

real*8, dimension(:), allocatable :: sedTrapDeployDepths
real*8, dimension(:), allocatable :: imagingDeployDepths

real*8, dimension(:,:,:), allocatable :: sedTrap
real*8, dimension(:,:,:), allocatable :: sedTrapSf  

real*8, dimension(:,:,:), allocatable, target :: particleFlux   ! mg m-2 d-1
real*8, dimension(:,:,:), allocatable, target :: particleFluxSf ! mg m-2 d-1

real*8, dimension(:), allocatable :: cameraUpperBoundaries, &
									 cameraLowerBoundaries, &
									 cameraSampledDepth
									 
real*8, dimension(:,:), allocatable :: cameraSampledVolume ! L
								
real*8, dimension(:,:), allocatable :: lossTerms	
 
real*8, dimension(:,:,:,:), allocatable :: attsInSizeClass
real*8, dimension(:,:,:,:), allocatable :: attsInVeloClass	
real*8, dimension(:,:,:), allocatable   :: attsInVolumeClass	
real*8, dimension(:,:,:,:), allocatable :: attsInSizeClassSf
real*8, dimension(:,:,:,:), allocatable :: attsInVeloClassSf
real*8, dimension(:,:,:,:), allocatable :: attsInSizeClassInst
real*8, dimension(:,:,:,:), allocatable :: attsInVeloClassInst	
real*8, dimension(:,:,:,:), allocatable :: attsInSizeClassSfInst
real*8, dimension(:,:,:,:), allocatable :: attsInVeloClassSfInst	

real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInSizeClass
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInVeloClass	
real*8, dimension(:,:,:), allocatable, target   :: avgAttsInVolumeClass
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInSizeClassSf
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInVeloClassSf
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInSizeClassInst
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInVeloClassInst
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInSizeClassSfInst
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInVeloClassSfInst

real*8, dimension(:,:,:,:), allocatable :: attsInMainType		
real*8, dimension(:,:,:,:), allocatable :: attsInMainTypeSf
real*8, dimension(:,:,:,:), allocatable :: attsInMainTypeInst
real*8, dimension(:,:,:,:), allocatable :: attsInMainTypeSfInst

real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInMainType
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInMainTypeSf
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInMainTypeInst
real*8, dimension(:,:,:,:), allocatable, target :: avgAttsInMainTypeSfInst

contains

! ========================================================================================	

subroutine ParticleCollectionVariables(nProfiles)

	integer, intent(in) :: nProfiles
	logical, dimension(:), allocatable :: matchCriteria
	integer :: i, iImagingDeployDepth, iProfile
	character(len=30) :: fmt

	fmt = '(7F9.1)'

#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS

	! Configuration under fixed initial particle conditions
	nSedTrapDeployDepths = 1
	nImagingDeployDepths = 1
	allocate(sedTrapDeployDepths(nSedTrapDeployDepths))
	allocate(imagingDeployDepths(nImagingDeployDepths))
	sedTrapDeployDepths(:) = 50d0 ! fixed depth for sediment traps
	imagingDeployDepths(:) = 50d0 ! fixed depth for imaging devices
	allocate(lossTerms(2,nVolumeClasses))	
	lossTerms(:,:) = 0d0
	
#else

	! ------ Sediment trap depths: match sediment trap depths within valid depth range ------
	allocate(matchCriteria(maxNumSedTrapDeployDepths))
	matchCriteria(:) = (possibleSedTrapDeployDepths(:) >= ztop(1)) &
		.and. (possibleSedTrapDeployDepths(:) < zbot(nDepthLayers)) 	
	nSedTrapDeployDepths = COUNT(matchCriteria(:))
	allocate(sedTrapDeployDepths(nSedTrapDeployDepths))
	sedTrapDeployDepths = (PACK(possibleSedTrapDeployDepths(:), matchCriteria(:)))
	deallocate(matchCriteria)
	
	print *
	print *, 'Sediment trap deployment depths (m):'
	print fmt, (sedTrapDeployDepths(i), i=1,nSedTrapDeployDepths)

	! ------ Imaging depths: match imaging device depths within valid depth range ------
	allocate(matchCriteria(maxNumImagingDeployDepths))
	matchCriteria(:) = (possibleImagingDeployDepths(:) >= ztop(1)) &
		.and. (possibleImagingDeployDepths(:) < zbot(nDepthLayers))	
	nImagingDeployDepths = COUNT(matchCriteria(:))
	allocate(imagingDeployDepths(nImagingDeployDepths))
	imagingDeployDepths(:) = (PACK(possibleImagingDeployDepths(:), matchCriteria(:)))
	deallocate(matchCriteria)
	
	print *
	print *, 'Imaging deployment depths (m):'
	print fmt, (imagingDeployDepths(i), i=1,nImagingDeployDepths)
	print *

#endif

	! ------ Allocate and compute imaging device boundary depths and sampled volumes ------
	allocate(cameraUpperBoundaries(nImagingDeployDepths))
	allocate(cameraLowerBoundaries(nImagingDeployDepths))
	allocate(cameraSampledDepth(nImagingDeployDepths))
	allocate(cameraSampledVolume(nImagingDeployDepths,nProfiles))

	cameraUpperBoundaries(1) = 0d0
	do iImagingDeployDepth = 1, nImagingDeployDepths
		! Handle the last imaging depth separately
		if (iImagingDeployDepth == nImagingDeployDepths) then		
			cameraLowerBoundaries(iImagingDeployDepth) = imagingDeployDepths(iImagingDeployDepth)
			cameraSampledDepth(iImagingDeployDepth) = cameraLowerBoundaries(iImagingDeployDepth) &
				- cameraUpperBoundaries(iImagingDeployDepth) ! m
			exit
		end if
		! For other depths, assign upper and lower boundaries
		cameraUpperBoundaries(iImagingDeployDepth+1) = imagingDeployDepths(iImagingDeployDepth)
		cameraLowerBoundaries(iImagingDeployDepth) = imagingDeployDepths(iImagingDeployDepth)
		cameraSampledDepth(iImagingDeployDepth) = cameraLowerBoundaries(iImagingDeployDepth) &
			- cameraUpperBoundaries(iImagingDeployDepth) ! m
	end do	
	
	! ------ Compute sampled volume for each profile ------
	do iProfile = 1, nProfiles
		cameraSampledVolume(:,iProfile) = cameraSampledDepth(:)*gridCellArea(iProfile)*1d3 ! L
	end do
	
	print *, 'Camera upper boundaries'
	print fmt, (cameraUpperBoundaries(i), i=1,nImagingDeployDepths)
	print *, 'Camera lower boundaries'
	print fmt, (cameraLowerBoundaries(i), i=1,nImagingDeployDepths)
	print *, 'Camera sampled depths (m):'
	print fmt, (cameraSampledDepth(i), i=1,nImagingDeployDepths)
	print *, 'Camera sampled volume (L):'
	print fmt, (cameraSampledVolume(i,1), i=1,nImagingDeployDepths)
	print *
    
	! ------ Allocate space for arrays for particle flux and attributes ------
	allocate(sedTrap(5,nSedTrapDeployDepths,nProfiles))
	allocate(sedTrapSf(5,1,nProfiles))
	
	allocate(particleFlux(5,nSedTrapDeployDepths,nProfiles))
	allocate(particleFluxSf(5,1,nProfiles))
	
	allocate(attsInSizeClass(nSizeClasses,nImagingDeployDepths,14,nProfiles))
	allocate(attsInVeloClass(nVeloClasses,nImagingDeployDepths,14,nProfiles))
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	allocate(attsInVolumeClass(nVolumeClasses,11,nProfiles))
#endif
	allocate(attsInSizeClassSf(nSizeClasses,1,14,nProfiles))
	allocate(attsInVeloClassSf(nVeloClasses,1,14,nProfiles))
	allocate(attsInSizeClassInst(nSizeClasses,nImagingDeployDepths,14,nProfiles))
	allocate(attsInVeloClassInst(nVeloClasses,nImagingDeployDepths,14,nProfiles))
	allocate(attsInSizeClassSfInst(nSizeClasses,1,14,nProfiles))
	allocate(attsInVeloClassSfInst(nVeloClasses,1,14,nProfiles))

	allocate(avgAttsInSizeClass(nSizeClasses,nImagingDeployDepths,14,nProfiles))
	allocate(avgAttsInVeloClass(nVeloClasses,nImagingDeployDepths,14,nProfiles))
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	allocate(avgAttsInVolumeClass(nVolumeClasses,11,nProfiles))
#endif
	allocate(avgAttsInSizeClassSf(nSizeClasses,1,14,nProfiles))
	allocate(avgAttsInVeloClassSf(nVeloClasses,1,14,nProfiles))
	allocate(avgAttsInSizeClassInst(nSizeClasses,nImagingDeployDepths,14,nProfiles))
	allocate(avgAttsInVeloClassInst(nVeloClasses,nImagingDeployDepths,14,nProfiles))
	allocate(avgAttsInSizeClassSfInst(nSizeClasses,1,14,nProfiles))
	allocate(avgAttsInVeloClassSfInst(nVeloClasses,1,14,nProfiles))

	allocate(attsInMainType(nMainParticleTypes,nImagingDeployDepths,7,nProfiles))
	allocate(attsInMainTypeSf(nMainParticleTypes,1,7,nProfiles))
	allocate(attsInMainTypeInst(nMainParticleTypes,nImagingDeployDepths,7,nProfiles))
	allocate(attsInMainTypeSfInst(nMainParticleTypes,1,7,nProfiles))
				
	allocate(avgAttsInMainType(nMainParticleTypes,nImagingDeployDepths,7,nProfiles))
	allocate(avgAttsInMainTypeSf(nMainParticleTypes,1,7,nProfiles))
	allocate(avgAttsInMainTypeInst(nMainParticleTypes,nImagingDeployDepths,7,nProfiles))
	allocate(avgAttsInMainTypeSfInst(nMainParticleTypes,1,7,nProfiles))

	! ------ Initialise particle flux and attributes arrays to zero	 ------
	sedTrap(:,:,:)   = 0d0
	sedTrapSf(:,:,:) = 0d0
					
	attsInSizeClass(:,:,:,:)       = 0d0
	attsInVeloClass(:,:,:,:)       = 0d0
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
	attsInVolumeClass(:,:,:)       = 0d0
#endif		
	attsInSizeClassSf(:,:,:,:)     = 0d0
	attsInVeloClassSf(:,:,:,:) 	   = 0d0
	attsInSizeClassInst(:,:,:,:)   = 0d0
	attsInVeloClassInst(:,:,:,:)   = 0d0	
	attsInSizeClassSfInst(:,:,:,:) = 0d0
	attsInVeloClassSfInst(:,:,:,:) = 0d0
	
	attsInMainType(:,:,:,:)        = 0d0
	attsInMainTypeSf(:,:,:,:)      = 0d0
	attsInMainTypeInst(:,:,:,:)    = 0d0
	attsInMainTypeSfInst(:,:,:,:)  = 0d0

end subroutine ParticleCollectionVariables

! ========================================================================================	

end module modelparticlecollection