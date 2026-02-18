module modelgrid

! ----------------------------------------------------------------------------------------
! This module is concerned with carrying the model physical grid.
! ----------------------------------------------------------------------------------------

use modelparameters, only: nDepthLayers, filenameZub, filenameZlb

implicit none
external :: read_r8_field
public

real*8, dimension(:), allocatable :: gridCellArea							   
real*8, dimension(:), allocatable :: ztop, zbot, zmid

contains

! ========================================================================================	

subroutine AreaGrid(nProfiles)

	integer, intent(in) :: nProfiles

	allocate(gridCellArea(nProfiles))

    gridCellArea(:) = 0.10d0 ! m2, initialise to a common value for all profiles
		
end subroutine AreaGrid

! ========================================================================================	

subroutine DepthLayerGrid( )

	integer :: i
	character(len=30) :: fmt

	fmt = '(7F9.1)'

	allocate(ztop(nDepthLayers))
	allocate(zbot(nDepthLayers))
	allocate(zmid(nDepthLayers))

! 	E.g.: 0, 10, 20, 30, 40, ..., 4480, 4490	
	call read_r8_field(nDepthLayers,1,1,ztop,filenameZub)
!	E.g.: 10, 20, 30, 40, 50, ..., 4490, 4500
	call read_r8_field(nDepthLayers,1,1,zbot,filenameZlb)
	
	! Compute mid-depth for each layer
	zmid = zbot - 0.5d0*(zbot - ztop)
	
end subroutine DepthLayerGrid

! ========================================================================================

end module modelgrid