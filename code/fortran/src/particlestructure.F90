module particlestructure

! ----------------------------------------------------------------------------------------
! This module is concerned with initialising the particle array.
! ----------------------------------------------------------------------------------------

use iso_c_binding, only: C_INT, C_DOUBLE

implicit none
public

type, bind(c) :: lagrangianStateVars ! array of structures
#include "pstructc.h"
end type lagrangianStateVars

contains

! ========================================================================================

subroutine InitialiseParticleArray(particle, nClusters)

	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	integer, intent(in) :: nClusters
	
	integer :: iCluster	

	do iCluster = 1, nClusters
		call InitialiseParticle(particle(iCluster))
	end do	

end subroutine InitialiseParticleArray

! ========================================================================================

subroutine InitialiseParticle(particle) bind(c, name='InitialiseParticle')

	type(lagrangianStateVars), intent(inout) :: particle
	
	particle%id            = 0
	particle%initType      = 0
	particle%initPft       = 0
	particle%phase         = 0
	particle%living        = 0
	particle%faecal        = 0
	particle%tstepCreat    = 0
	
	particle%nPxC  		   = 0d0
	particle%nPpxP 		   = 0d0
	particle%nPpxC 		   = 0d0

	particle%molesOrgC     = 0d0
	particle%molesTepC     = 0d0
	particle%molesMineral  = 0d0		
	particle%massOrgMatter = 0d0
	particle%massTep       = 0d0
	particle%massMineral   = 0d0
	particle%depth         = 0d0
	particle%mass          = 0d0
	particle%stickiness    = 0d0
	particle%fracDim       = 0d0
	particle%radiusPp      = 0d0
	particle%radius        = 0d0
	particle%solidVolume   = 0d0
	particle%porosity      = 0d0
	particle%density       = 0d0
	particle%excessDensity = 0d0
	particle%velocity      = 0d0
	
end subroutine InitialiseParticle

! ========================================================================================
	
end module particlestructure