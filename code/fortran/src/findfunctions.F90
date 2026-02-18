module findfunctions

! ----------------------------------------------------------------------------------------
! This module contains functions that find objects.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use modelparameters, only: nVolumeClasses, particleVolumeClasses
use sanitychecks, only: WriteStatusAndStop

implicit none
private
public :: FindActiveClusters, FindPhase1clusters, FindClustersInDepthLayer, &
	FindDepthLayerMidpointDepth, FindDepthLayerIndex, FindIndexToNearestLayer, &
	FindParticleAttributeClass, FindParticleVolumeClass

contains

! ========================================================================================

subroutine FindActiveClusters(nActiveClusters, activeClusterIndices, particle, nClusters, iLastLocus)

	!-------------------------------------------------------------------------------------
	! These are clusters with an ID, and can be phase = 1, 2 or 3
	!-------------------------------------------------------------------------------------	
	
	integer, intent(in) :: nClusters, iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	integer, intent(out) :: nActiveClusters
	integer, dimension(:), allocatable, intent(out) :: activeClusterIndices
	
	integer, dimension(:), allocatable :: indexArray
	logical, dimension(:), allocatable :: maskMatch
	integer :: i

	! ------ Sanity checks first ------
	if (iLastLocus < 0 .or. iLastLocus > nClusters) then
  		write(*,*) 'ERROR: in FindActiveClusters we found an invalid iLastLocus', iLastLocus, 'nClusters=', nClusters
  		call WriteStatusAndStop( )
	end if
	if (any(.not. ieee_is_finite(dble(particle(1:iLastLocus)%id)))) then
    	write(*,*) 'ERROR: non-finite particle IDs in FindActiveClusters'
    	call WriteStatusAndStop()
	end if
	
	allocate(maskMatch(iLastLocus), indexArray(iLastLocus))
                                                    ! e.g. for iLastLocus = 10, particle(1:iLastLocus)%id = (1, 0, 99, 5, 5, 0, 20, 10, 9, 9)
	maskMatch(:) = particle(1:iLastLocus)%id > 0  !                                               maskMatch = (1, 0, 1,  1, 1, 0,  1,  1, 1, 1)
	indexArray(:) = [(i, i = 1, iLastLocus, 1)]     !                                            indexArray = (1, 2, 3,  4, 5, 6,  7,  8, 9, 10)

	nActiveClusters = COUNT(maskMatch)             !                                          nActiveClusters = 8
	if (allocated(activeClusterIndices)) deallocate(activeClusterIndices)
	allocate(activeClusterIndices(nActiveClusters))

	activeClusterIndices(:) = PACK(indexArray(:), maskMatch(:)) !                        activeClusterIndices = (1, 3, 4, 5, 7, 8, 9, 10)

	deallocate(maskMatch, indexArray)

end subroutine FindActiveClusters

! ========================================================================================

subroutine FindPhase1clusters(nPhase1clusters, phase1clusterIndices, particle, nClusters, iLastLocus)

	!-------------------------------------------------------------------------------------
	! These are clusters with an ID and phase = 1 only
	!-------------------------------------------------------------------------------------	
	
	integer, intent(in) :: nClusters, iLastLocus
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	integer, intent(out) :: nPhase1clusters
	integer, dimension(:), allocatable, intent(out) :: phase1clusterIndices
	
	integer, dimension(:), allocatable :: indexArray
	logical, dimension(:), allocatable :: maskMatch
	integer :: i
	
	! ------ Sanity checks first ------
	if (iLastLocus < 0 .or. iLastLocus > nClusters) then
  		write(*,*) 'ERROR: in FindPhase1clusters we found an invalid iLastLocus', iLastLocus, 'nClusters=', nClusters
  		call WriteStatusAndStop( )
	end if
	if (any(.not. ieee_is_finite(dble(particle(1:iLastLocus)%id)))) then
    	write(*,*) 'ERROR: non-finite particle IDs in FindPhase1clusters'
    	call WriteStatusAndStop()
	end if
	if (any(particle(1:iLastLocus)%phase < 0 .or. particle(1:iLastLocus)%phase > 3)) then
    	write(*,*) 'ERROR: invalid particle phase detected in FindPhase1clusters'
    	call WriteStatusAndStop()
	end if

	allocate(maskMatch(iLastLocus), indexArray(iLastLocus))

	maskMatch(:) = particle(1:iLastLocus)%phase == 1 ! e.g. for iLastLocus = 10, maskMatch = (1, 0, 1, 0, 1, 0, 0, 0, 1, 0)
	indexArray(:) = [(i, i = 1, iLastLocus, 1)]      !                          indexArray = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
	
	nPhase1clusters = COUNT(maskMatch)               !                     nPhase1clusters = 4
	if (allocated(phase1clusterIndices)) deallocate(phase1clusterIndices)
	allocate(phase1clusterIndices(nPhase1clusters))
  
	phase1clusterIndices(:) = PACK(indexArray(:), maskMatch(:)) !     phase1clusterIndices = (1, 3, 5, 9)

	deallocate(maskMatch, indexArray)

end subroutine FindPhase1clusters

! ========================================================================================

subroutine FindClustersInDepthLayer(nLayerClusters, layerClusterIndices, particle, &
	nClusters, nActiveClusters, activeClusterIndices, ztopLayer, zbotLayer)

	integer, intent(in) :: nClusters, nActiveClusters
	integer, dimension(nActiveClusters), intent(in) :: activeClusterIndices	
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	real*8, intent(in) :: ztopLayer, zbotLayer
	integer, intent(out) :: nLayerClusters
	integer, dimension(:), allocatable, intent(out) :: layerClusterIndices	

	logical, dimension(:), allocatable :: maskMatch
	
	! ------ Sanity checks first ------
	if (.not. ieee_is_finite(ztopLayer) .or. .not. ieee_is_finite(zbotLayer) .or. zbotLayer <= ztopLayer) then
    	write(*,*) 'ERROR: invalid depth bounds in FindClustersInDepthLayer'
    	call WriteStatusAndStop()
	end if
	if (any(activeClusterIndices < 1) .or. any(activeClusterIndices > nClusters)) then
    	write(*,*) 'ERROR: invalid activeClusterIndices in FindClustersInDepthLayer'
    	call WriteStatusAndStop()
	end if
			
	allocate(maskMatch(nActiveClusters))
											                      
	maskMatch(:) = particle(activeClusterIndices)%phase == 1 &          ! activeClusterIndices = (1, 3, 4, 5, 7, 8, 9, 10)
		     .and. particle(activeClusterIndices)%depth >= ztopLayer &  !      nActiveClusters = 8
		     .and. particle(activeClusterIndices)%depth < zbotLayer     !            maskMatch = (1, 0, 0, 1, 0, 1, 0, 0)

	nLayerClusters = COUNT(maskMatch)                                   !       nLayerClusters = 3
	if (allocated(layerClusterIndices)) deallocate(layerClusterIndices)
	allocate(layerClusterIndices(nLayerClusters))
	  
	layerClusterIndices(:) = PACK(activeClusterIndices(:), maskMatch(:))!  layerClusterIndices = (1, 5, 8)

	deallocate(maskMatch)
		
end subroutine FindClustersInDepthLayer
		
! ========================================================================================

function FindDepthLayerMidpointDepth(iDepth, nLayers, depthUpper, depthLower) result(zmid)

	real*8 :: zmid
	integer, intent(in) :: nLayers, iDepth
	real*8, dimension(nLayers), intent(in) :: depthUpper, depthLower
	
	! ------ Sanity checks first ------
	if (depthLower(iDepth) < depthUpper(iDepth)) then
    	write(*,*) 'ERROR: inverted layer bounds in FindDepthLayerMidpointDepth'
    	call WriteStatusAndStop()
	end if
 	if (iDepth < 1 .or. iDepth > nLayers) then
  		write(*,*) 'ERROR: in FindDepthLayerMidpointDepth, we found an iDepth out-of-range', iDepth
  		call WriteStatusAndStop( )
	end if
	
	zmid = depthLower(iDepth) - 0.5d0*(depthLower(iDepth)-depthUpper(iDepth))
	
end function FindDepthLayerMidpointDepth

! ========================================================================================

function FindDepthLayerIndex(particleDepth, nLayers, depthUpper, depthLower) result(iDepthLayer)
    
    integer :: iDepthLayer
    integer, intent(in) :: nLayers
    real*8, intent(in) :: particleDepth
    real*8, dimension(nLayers), intent(in) :: depthUpper, depthLower
    integer :: i

	! ------ Sanity checks first ------
    if (.not. ieee_is_finite(particleDepth) .or. particleDepth < 0d0) then 
    	write(*,*) 'ERROR: invalid particleDepth detected when finding its depth layer', particleDepth
        call WriteStatusAndStop( )
    end if
    if (nLayers < 1) then
    	write(*,*) 'ERROR: nLayers < 1 in FindDepthLayerIndex'
    	call WriteStatusAndStop()
	end if
	if (any(depthLower < depthUpper)) then
    	write(*,*) 'ERROR: invalid layer geometry in FindDepthLayerIndex'
    	call WriteStatusAndStop()
	end if
    
    iDepthLayer = 0
    do i = 1, nLayers
    	if (i < nLayers) then
      		! normal layers: include upper bound, exclude lower bound
      		if (particleDepth >= depthUpper(i) .and. particleDepth < depthLower(i)) then
        		iDepthLayer = i
        		return
      		end if
    	else
      		! last layer: make lower boundary inclusive so bottom-most point maps in
      		if (particleDepth >= depthUpper(i) .and. particleDepth <= depthLower(i)) then
        		iDepthLayer = i
        		return
      		end if
    	end if
  	end do
  	
  	! Boundary fallback (rare): map very shallow or very deep values to ends
  	if (particleDepth < depthUpper(1)) then
    	iDepthLayer = 1
  	else if (particleDepth > depthLower(nLayers)) then
    	iDepthLayer = nLayers
  	end if
  	if (iDepthLayer < 1 .or. iDepthLayer > nLayers) then
    	write(*,*) 'ERROR: failed to assign depth layer for depth', particleDepth
    	call WriteStatusAndStop()
	end if
    
end function FindDepthLayerIndex

! ========================================================================================

function FindIndexToNearestLayer(targetDepth,nLocalDepthLayers, zmid) result(idx)

	real*8, intent(in) :: targetDepth
	integer, intent(in) :: nLocalDepthLayers
	real*8, dimension(nLocalDepthLayers), intent(in) :: zmid
	
	integer :: idx
	real*8 :: d, dmin
	integer :: k

	idx  = 1
	dmin = ABS(zmid(1) - targetDepth)
	do k = 2, nLocalDepthLayers
		d = ABS(zmid(k) - targetDepth)
		if (d < dmin) then
			dmin = d
			idx  = k
		end if
	end do

end function FindIndexToNearestLayer

! ========================================================================================

function FindParticleAttributeClass(nAttributeClasses, listOfAttributeClasses, particleAttribute) result(classIdx)

	!-------------------------------------------------------------------------------------
	! Classify the particle in a diameter or velocity class according to already predefined 
	! categories.
	!-------------------------------------------------------------------------------------
	
	integer :: classIdx
	integer, intent(in) :: nAttributeClasses
	real*8, intent(in) :: particleAttribute
	real*8, dimension(nAttributeClasses), intent(in) :: listOfAttributeClasses

	integer :: i
  	
    if (particleAttribute < listOfAttributeClasses(1)) then
        classIdx = 1
    else if (particleAttribute > listOfAttributeClasses(nAttributeClasses)) then
        classIdx = nAttributeClasses
	else
		classIdx = nAttributeClasses
		do i = 1, nAttributeClasses
        	if (particleAttribute <= listOfAttributeClasses(i)) then
            	classIdx = i
            	exit
        	end if
   	 	end do
	end if

end function FindParticleAttributeClass

! ========================================================================================

function FindParticleVolumeClass(particleVolume) result(volClassIdx)

	!-------------------------------------------------------------------------------------
	! Classify particles into volume sections as in Jackson's model.
	!-------------------------------------------------------------------------------------
	
	integer :: volClassIdx
	real*8, intent(in) :: particleVolume
	
	integer :: i

	volClassIdx = 0
	do i = 1, nVolumeClasses
		if (particleVolume <= particleVolumeClasses(i)) then
			volClassIdx = i
			exit
		end if	
	end do
	
	if (volClassIdx == 0) then
    	write(*,*) 'ERROR: Particle volume class not assigned correctly.'
		call WriteStatusAndStop( )
	end if

end function FindParticleVolumeClass

! ========================================================================================

end module findfunctions