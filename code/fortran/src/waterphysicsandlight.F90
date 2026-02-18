module waterphysicsandlight

! ----------------------------------------------------------------------------------------
! This module handles functions that calculate physical and optical properties in seawater.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use modelconstants, only: PI
use sanitychecks, only: WriteStatusAndStop
use montecarlosampling, only: RandProbCase

implicit none
private
public ::  WaterKinematicViscosity, TurbulentKineticEnergyDissipationRateDepthProfile, ParProfile

contains

! ========================================================================================

function WaterKinematicViscosity(dynvisco, dens)

	real*8 :: WaterKinematicViscosity ! denoted by the Greek letter nu
	real*8, intent(in) :: dynvisco, dens

	WaterKinematicViscosity = dynvisco/dens ! cm2 s-1
	WaterKinematicViscosity = WaterKinematicViscosity * 1d-4 ! m2 s-1

end function WaterKinematicViscosity

! ========================================================================================

subroutine TurbulentKineticEnergyDissipationRateDepthProfile(tkeDepthProfile, &
	nLocalDepthLayers, localSeafloorDepth)

	integer, intent(in) :: nLocalDepthLayers
	real*8, intent(in) :: localSeafloorDepth
	real*8, dimension(nLocalDepthLayers), intent(out) :: tkeDepthProfile
	
	real*8 :: minTkeRange, maxTkeRange, mu, sigma, tkeSurface
	real*8, dimension(:), allocatable :: tkeRange, logNormalPdf
	integer :: n, n1, i, iProbTke

	! TKE ranges from Kiorboe & Saiz (1995), surface values
	
	if (localSeafloorDepth <= 200d0) then ! shelf seas and coastal ocean
		minTkeRange = -7 ! 1d-7 m2 s-3 (only write the power)
		maxTkeRange = -4 ! 1d-4 m2 s-3 (only write the power)
	else ! open ocean
		minTkeRange = -9 ! 1d-8 m2 s-3 (only write the power)
		maxTkeRange = -6  ! 1d-6 m2 s-3 (only write the power)	
	end if
	n = FLOOR(10**(ABS(minTkeRange) - ABS(maxTkeRange)))
	allocate(tkeRange(n), logNormalPdf(n))

	! Construct a log-normal distribution curve with the values from the chosen range of TKE
	! and randomly sample it

	! (1) logspace of the range of tke values (edit logspace)
	n1 = n-1
	tkeRange(:) = (/( (minTkeRange + i*((maxTkeRange-minTkeRange)/n1)), i=0, n1 )/) 
	tkeRange(:) = 10d0**tkeRange(:)
	
	! (2) mean and variance
	mu = SUM(LOG(tkeRange(:)))/n  ! ln of the mean value in the tke range
 	sigma = 1.2d0
 	
	! (3) log-normal probability distribution function (PDF)
	logNormalPdf(:) = EXP( -0.5d0 * ((LOG(tkeRange(:))-mu)/sigma)**2 ) / (tkeRange(:)*sigma*SQRT(2d0*PI))
	if (.not. ieee_is_finite(sum(logNormalPdf)) .or. sum(logNormalPdf) <= 0d0) then
    	write(*,*) 'ERROR: invalid TKE PDF normalisation'
    	call WriteStatusAndStop()
	end if
	logNormalPdf(:) = logNormalPdf(:) / SUM(logNormalPdf(:))
	if (abs(sum(logNormalPdf) - 1d0) > 1d-12) then
    	write(*,*) 'ERROR: TKE PDF does not sum to 1'
    	call WriteStatusAndStop()
	end if
	!print *, 'Is SUM(logNormalPdf(:))=1?', SUM(logNormalPdf(:))
	
	! (4) using the cumulative distribution function (CDF), randomly choose a TKE value 
	iProbTke = RandProbCase(n, logNormalPdf(:))
	if (iProbTke < 1 .or. iProbTke > n) then
    	write(*,*) 'ERROR: RandProbCase returned out-of-bounds index for iProbTke'
    	call WriteStatusAndStop()
	end if
	tkeSurface = tkeRange(iProbTke) ! m2 s-3
	if (.not. ieee_is_finite(tkeSurface) .or. tkeSurface <= 0d0) then
    	write(*,*) 'ERROR: invalid sampled TKE value'
    	call WriteStatusAndStop()
	end if

	! The vertical distribution of TKE is constant
	tkeDepthProfile(:) = tkeSurface
	
	deallocate(tkeRange)
	deallocate(logNormalPdf)
	
end subroutine TurbulentKineticEnergyDissipationRateDepthProfile

! ========================================================================================

subroutine ParProfile(parz, nLocalDepthLayers, z, parO, chla, mld)

	!-------------------------------------------------------------------------------------
	! Equations as in Fox et al. (2024), after Morel & Maritorena (2001) and Morel et al. (2007)
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nLocalDepthLayers
	real*8, intent(in) :: parO, chla, mld
	real*8, intent(in), dimension(nLocalDepthLayers) :: z
	real*8, dimension(nLocalDepthLayers), intent(out) :: parz
	
	integer :: i
    real*8 :: kd, kpar

	! Calculate diffuse attenuation coefficient at 490 nm (kd)
    kd = 0.0166d0 + 0.07242d0 * (chla**0.68955d0)

	! Calculate diffuse attenuation coefficient for Photosynthetically Available Radiation (kPAR)
	! based on Mixed Layer Depth (MLD)
    if (mld <= 1.0d0 / kd) then
        kpar = 0.0864d0 + 0.884d0 * kd - (0.00137d0 / kd)
    else
        kpar = 0.0665d0 + 0.874d0 * kd - (0.00121d0 / kd)
    end if

	! Calculate PAR(z) for each depth in the array
    do i = 1, nLocalDepthLayers
        parz(i) = parO * EXP(-kpar * z(i)) ! W m-2
    end do

	!kpar = SCATTERING_COEFF_SEAWATER + (ABSORPTION_COEFF_CHLA*chla)
	!ParProfile = parO * EXP( -kpar*z ) ! W m-2

end subroutine ParProfile

! ========================================================================================

end module waterphysicsandlight