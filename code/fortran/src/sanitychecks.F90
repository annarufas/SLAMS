module sanitychecks

! ----------------------------------------------------------------------------------------
! This module is concerned with sanity checks.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use modelparameters, only: detection_limit_poc, detection_limit_calc, detection_limit_opal, &
	detection_limit_clay, iCalcite, iOpal, iClay

implicit none
public :: IsQuantityEffectivelyZero, NearlyEqual, IsExceedingInitialAmount, &
	CheckParticleSanity, PrintParticleProperties, WriteStatus, WriteStatusAndStop

contains

! ========================================================================================

logical function IsQuantityEffectivelyZero(x, analyticalDetectionLimit)

    !-------------------------------------------------------------------------------------
  	! Returns .true. if |x| <= analyticalDetectionLimit. Good for decisions of the type:
  	! "is a quantity small enough to treat as zero for the model"?
  	!-------------------------------------------------------------------------------------

    real*8, intent(in) :: x, analyticalDetectionLimit

    if (.not. ieee_is_finite(x)) then
		IsQuantityEffectivelyZero = .false.
      	return
    end if

    IsQuantityEffectivelyZero = (abs(x) <= analyticalDetectionLimit)
    
end function IsQuantityEffectivelyZero
  
! ========================================================================================

logical function NearlyEqual(a, b, analyticalDetectionLimit)

  	!-------------------------------------------------------------------------------------
  	! Returns .true. if a and b are equal detection limit. Good for decisions of the type:
  	! "are two quantities equal within detection limit"?
  	!-------------------------------------------------------------------------------------
  
    real*8, intent(in) :: a, b, analyticalDetectionLimit
    real*8 :: diff
    
     if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
      	NearlyEqual = .false.
      	return
    end if
    
    diff = a - b
    NearlyEqual = IsQuantityEffectivelyZero(diff, analyticalDetectionLimit)
    
end function NearlyEqual

! ========================================================================================

logical function IsExceedingInitialAmount(newValue, initialValue, analyticalDetectionLimit)

    real*8, intent(in) :: newValue, initialValue, analyticalDetectionLimit
    real*8 :: diff
   
    diff = newValue - initialValue
    IsExceedingInitialAmount = (diff > analyticalDetectionLimit) ! if diff is positive and greater than tol --> violation

end function IsExceedingInitialAmount
  
! ========================================================================================

subroutine CheckParticleSanity(particle, nClusters, iCluster, processDescription)

	!-------------------------------------------------------------------------------------
	! Returns .true. if particle passes sanity checks, .false. otherwise. Those checks
	! test whether the particle is physically consistent or not by looking at the mass
	! composition and the geometry, like the amount of primary particles. If the particle fails
	! the sanity check, the model stops running and creates a file that stops the core module.
	! This should not happen more than once.
	! The functionality is identical. I should use the write statements instead,
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	character(len=*), intent(in), optional :: processDescription
	
	logical :: ok
	character(len=:), allocatable :: procLabel
	
	if (present(processDescription)) then
      	procLabel = trim(processDescription)
    else
      	procLabel = 'unspecified'
    end if
	
	ok = .true.
	
	! 1) Particle counts must be >= 1
	if (particle(iCluster)%nPpxP < 1d0 .or. particle(iCluster)%nPxC < 1d0) then
	  	ok = .false.
	  	write(*,*) 'SANITY FAIL: nPpxP < 1 or nPxC < 1'
	  	write(*,*) '  nPpxP, nPxC, nPpxC =', particle(iCluster)%nPpxP, particle(iCluster)%nPxC, particle(iCluster)%nPpxC
	  	write(*,*) '  phase =', particle(iCluster)%phase
	  	write(*,*) '  mass OM=', particle(iCluster)%massOrgMatter
	end if
	
	! 2) Negative material: allow tiny negative numerical noise up to tolerance, but fail if < -tol
	if (particle(iCluster)%molesOrgC < -detection_limit_poc .or. &
		particle(iCluster)%molesTepC < -detection_limit_poc .or. &
		particle(iCluster)%molesMineral(iCalcite) < -detection_limit_calc .or. &
		particle(iCluster)%molesMineral(iOpal) < -detection_limit_opal .or. &
		particle(iCluster)%molesMineral(iClay) < -detection_limit_clay) then
	  	ok = .false.
	  	write(*,*) 'SANITY FAIL: negative material beyond numerical tolerance'
	  	write(*,*) '  molesOrgC, molesTepC =', particle(iCluster)%molesOrgC, particle(iCluster)%molesTepC
	  	write(*,*) '  molesMineral =', particle(iCluster)%molesMineral(:)
	end if
	
	! 3) NaN checks for important derived fields
	if (.not. ieee_is_finite(particle(iCluster)%massOrgMatter) .or. .not. ieee_is_finite(particle(iCluster)%massTep)) then
		ok = .false.
	  	write(*,*) 'SANITY FAIL: NaN mass fields'
	  	write(*,*) '  molesOrgC, molesTepC =', particle(iCluster)%molesOrgC, particle(iCluster)%molesTepC
	  	write(*,*) '  molesMineral =', particle(iCluster)%molesMineral(:)
	end if
	if (.not. ieee_is_finite(particle(iCluster)%velocity)) then
		ok = .false.
	  	write(*,*) 'SANITY FAIL: NaN velocity field'
	  	write(*,*) '  velocity =', particle(iCluster)%velocity
	  	write(*,*) '  molesOrgC, molesTepC =', particle(iCluster)%molesOrgC, particle(iCluster)%molesTepC
	  	write(*,*) '  molesMineral =', particle(iCluster)%molesMineral(:)
	end if
	
	! 4) Geometry check (rpp > rp)
	if (particle(iCluster)%radiusPp > particle(iCluster)%radius) then
	  	ok = .false.
	  	write(*,*) 'SANITY FAIL: radiusPp > radius'
	  	write(*,*) '  rpp, rp =', particle(iCluster)%radiusPp, particle(iCluster)%radius, 'Npp=', particle(iCluster)%nPpxP
	end if
	
	! If failed, write program_status.tmp and stop
	if (.not. ok) then
		write(*,*) 'FATAL: sanity failure for particle index', iCluster
		write(*,*) 'Process: ' // trim(procLabel)
		call WriteStatusAndStop( )
	end if
	
	! Maybe consider particle(iCluster)%phase = 3 for some physically inconsistent attributes?
	! Are there any circumstances when those might actually make sense?
	
end subroutine CheckParticleSanity

! ========================================================================================

subroutine PrintParticleProperties(particle, nClusters, iCluster, particleDescription)

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	character(len=*), intent(in), optional :: particleDescription
	
	character(len=:), allocatable :: partLabel
	
	if (present(particleDescription)) then
      	partLabel = trim(particleDescription)
    else
      	partLabel = 'unspecified'
    end if

	write(*,*)
	write(*,*) 'PROCESS: ' // trim(partLabel) 
	write(*,*) '  Phase', particle(iCluster)%phase
	write(*,*) '  Type', particle(iCluster)%initType
	write(*,*) '  Living', particle(iCluster)%living
	write(*,*) '  Faecal', particle(iCluster)%faecal
	write(*,*) '  Depth', particle(iCluster)%depth
	write(*,*) '  nPpxP', particle(iCluster)%nPpxP
	write(*,*) '  nPxC:', particle(iCluster)%nPxC
	write(*,*) '  Mass:', particle(iCluster)%mass
	write(*,*) '  Mol orgC (pmol/part):', particle(iCluster)%molesOrgC/1d-12
 	write(*,*) '  Mol TEPC (pmol/part):', particle(iCluster)%molesTepC/1d-12
	write(*,*) '  Mol mineral (pmol/part):', particle(iCluster)%molesMineral(:)/1d-12	
	write(*,*) '  Diameter', 2d0*particle(iCluster)%radius
	write(*,*) '  Porosity:', particle(iCluster)%porosity
	write(*,*) '  Stick:', particle(iCluster)%stickiness
	write(*,*) '  Density:', particle(iCluster)%density
	write(*,*) '  Velocity:', particle(iCluster)%velocity
	write(*,*) '  Vpor (um3):', particle(iCluster)%solidVolume/(1d0-particle(iCluster)%porosity)
	write(*,*) '  Vsol (um3):', particle(iCluster)%solidVolume
	write(*,*)

end subroutine PrintParticleProperties		

! ========================================================================================

subroutine WriteStatus(status)

    character(len=*), intent(in) :: status
    integer :: iu, ios

    open(newunit=iu, file='program_status.tmp', status='replace', action='write', iostat=ios)
    if (ios == 0) then
        write(iu,'(A)') trim(status)
        close(iu)
    else
        write(*,*) 'WARNING: could not open program_status.tmp (iostat=', ios, ')'
    end if
    
end subroutine WriteStatus

! ========================================================================================

subroutine WriteStatusAndStop(status)

    character(len=*), intent(in), optional :: status
    
    call WriteStatus(merge(status,'error',present(status)))
    stop 1
    
end subroutine WriteStatusAndStop

! ========================================================================================

end module sanitychecks