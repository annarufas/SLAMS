#include "blockdefinitions.h"

module calcparticleattributes

! ----------------------------------------------------------------------------------------
! This module handles the subroutines that update the particle attributes array.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars, InitialiseParticle
use safemath, only: SafeDivide, SafeFloorNonNegative
use modelconstants, only: MOLAR_MASS_CARBON, MOLAR_MASS_OPAL, MOLAR_MASS_CACO3, MOLAR_MASS_CLAY, &
	RHO_ORGMATTER, RHO_OPAL, RHO_CALCITE, RHO_CLAY, RHO_TEP, GRAVITY_CNT, REYNOLDS_LAMINAR_LIMIT, &
	SECONDS_PER_DAY, PI
use modelparameters, only: fractal_dimension_agg_max, fractal_dimension_agg_min, &
	C_frac_in_OM, C_frac_in_TEP, iOpal, iCalcite, iClay, Si2C_diat, Calc2C_cocco_max, &
	clay_quota, k_omega, choiceFractalDimensionScheme, choiceStickinessFunc, k_TEP, &
	choiceMineralEffOnFaecalPellPorosity, porosity_faecalpell_max, porosity_faecalpell_min, &
	stickiness_phyto_cell_initial, nPfts, operational_size_poc_min, detection_limit_poc, &
	detection_limit_calc, detection_limit_opal, detection_limit_clay
use sanitychecks, only: PrintParticleProperties, WriteStatusAndStop, CheckParticleSanity
	
implicit none
private
public :: ParticleDryMass, ParticleMaterialVolume, ParticleFractalDimension, ParticleRadius, &
	ParticlePorosity, FaecalPelletPorosity, ParticleDensity, ParticleStickiness, ParticleSettlingVelocity, &
	ComputeParticleAttributes, CompactParticleArray, ShrinkParticle, VolumetricFractionOfOrganicMatter, &
	ParticleUpdatedType, CalculateParticleReynoldsNumber, CalculateTurbulenceScaleReynoldsNumber
	
contains

! ========================================================================================

subroutine ParticleDryMass(particle, nClusters, iCluster)

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	real*8 :: massOrgC, massTepC

! 	Some definitions
! 	- Particulate organic carbon, C, is proportional (not equal) to organic matter, CHONP
! 	- Transparent exopolymer carbon, TEP-C, is proportional (not equal) to TEP matter.
!	    TEP matter is mainly carbohydrates (CHO) + a small fraction of other molecules 
!		+ lots of water (> 90% of TEP's weight)
!  	- 1 mol of Si = 1 mol of opal (SiO2·2(H2O))
!	- 1 mol of C (as C in calcite) = 1 mol of CaCO3	
!
!	Terminology 
! 	- Silicon: Si
! 	- Silica: SiO2
! 	- Silicate: SiO4
! 	- Silicic acid (silicate in dissolution): H4SiO4 (SiO2 + 2H2O)
! 	- Opal = biogenic silica = bSi (hydrated form of silica): SiO2·n(H2O)

	! Mass associated to organic carbon compounds    	
	massOrgC = particle(iCluster)%molesOrgC * MOLAR_MASS_CARBON ! g organic C per particle
	particle(iCluster)%massOrgMatter = massOrgC / C_frac_in_OM ! g organic matter per particle

	massTepC = particle(iCluster)%molesTepC * MOLAR_MASS_CARBON ! g C in a TEP aqueous particle
	particle(iCluster)%massTep = massTepC / C_frac_in_TEP ! g of TEP

	! Mass associated to inorganic compounds
	particle(iCluster)%massMineral(iOpal) = particle(iCluster)%molesMineral(iOpal) &
		* MOLAR_MASS_OPAL ! g bSi per particle
	particle(iCluster)%massMineral(iCalcite) = particle(iCluster)%molesMineral(iCalcite) &
		* MOLAR_MASS_CACO3 ! g CaCO3 per particle
	particle(iCluster)%massMineral(iClay) = particle(iCluster)%molesMineral(iClay) &
		* MOLAR_MASS_CLAY ! g clay per particle	

	! Particle dry mass		
	particle(iCluster)%mass = particle(iCluster)%massOrgMatter    &
							+ particle(iCluster)%massTep          &
							+ SUM(particle(iCluster)%massMineral(:)) ! g

	! Sanity check
	if (particle(iCluster)%molesOrgC < -detection_limit_poc .or. &
		particle(iCluster)%molesTepC < -detection_limit_poc .or. &
		particle(iCluster)%molesMineral(iCalcite) < -detection_limit_calc .or. &
		particle(iCluster)%molesMineral(iOpal) < -detection_limit_opal .or. &
		particle(iCluster)%molesMineral(iClay) < -detection_limit_clay) then
	  	write(*,*) 'SANITY FAIL: negative material beyond numerical tolerance'
	  	write(*,*) '  molesOrgC, molesTepC =', particle(iCluster)%molesOrgC, particle(iCluster)%molesTepC
	  	write(*,*) '  molesMineral =', particle(iCluster)%molesMineral(:)
	  	call WriteStatusAndStop( )
	end if
	
end subroutine ParticleDryMass

! ========================================================================================

subroutine ParticleMaterialVolume(particle, nClusters, iCluster)

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle

	particle(iCluster)%solidVolume = 1d12*( particle(iCluster)%massOrgMatter/RHO_ORGMATTER           &
										  + particle(iCluster)%massTep/RHO_TEP						 &
										  + particle(iCluster)%massMineral(iOpal)/RHO_OPAL           &				
										  + particle(iCluster)%massMineral(iCalcite)/RHO_CALCITE     &
										  + particle(iCluster)%massMineral(iClay)/RHO_CLAY ) ! cm3 --> um3

end subroutine ParticleMaterialVolume

! ========================================================================================

subroutine ParticleFractalDimension(particle, nClusters, iCluster, particleType)

	integer, intent(in) :: nClusters, iCluster, particleType
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	real*8 :: harvest
	
	! Sanity check on particleType
    if (particleType < 1 .or. particleType > 9) then ! particleType in SLAMS-20 goes 1..9
    	write(*,*) 'SANITY FAIL: invalid particleType detected when calculating fractal dimension.'
        call WriteStatusAndStop( )
    end if
	
	! For aggregates (faecal and non-faecal)
	if (particleType <= 2) then
		select case (choiceFractalDimensionScheme)	
		case (1) ! fixed value
			particle(iCluster)%fracDim = (fractal_dimension_agg_max+fractal_dimension_agg_min)*0.5d0
		case (2) ! random value	
			call random_number( harvest )
			particle(iCluster)%fracDim = (fractal_dimension_agg_max-fractal_dimension_agg_min)*harvest + fractal_dimension_agg_min
		end select
		
	! For perfectly rounded particles (living phyto, dead phyto, TEP, clay, dead zoo, inorganic non-faecal single particles, organic non-faecal single particles)	
	else if (particleType > 2) then
		particle(iCluster)%fracDim = 3d0

	end if
	
end subroutine ParticleFractalDimension

! ========================================================================================

subroutine ParticleRadius(particle, nClusters, iCluster)

	!-------------------------------------------------------------------------------------
	! This subroutine allows rpp to exceed the radius normally used to classify a primary
	! particle as "particulate" material. This acts as a free boundary condition, avoiding
	! the need to artificially modify the number of primary particles in an aggregate in
	! order to preserve rpp at values ≥ 0.1 µm (or 0.2 ESD, the min. operational size of POC).
	! Primary particles that detach (e.g., due to fragmentation) and fall below the
	! detection threshold of 0.1 µm are automatically treated as dissolved material.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle

	real*8 :: primParticleMatVol
	
	! Sanity check before
	if (.not. ieee_is_finite(particle(iCluster)%solidVolume) .or. particle(iCluster)%solidVolume == 0d0) then
    	write(*,*) 'ERROR: invalid solidVolume in ParticleRadius', particle(iCluster)%solidVolume
    	write(*,*) '  mass:', particle(iCluster)%mass
    	write(*,*) '  nPxC, nPpxP:', particle(iCluster)%nPxC, particle(iCluster)%nPpxP
    	write(*,*) '  phase:', particle(iCluster)%phase
    	call WriteStatusAndStop( )
  	end if
  	
  	primParticleMatVol = SafeDivide(particle(iCluster)%solidVolume, particle(iCluster)%nPpxP, 0d0) ! um3
	if (primParticleMatVol <= 0d0) then
        write(*,*) 'ERROR: non-positive primParticleMatVol in ParticleRadius'
        call WriteStatusAndStop()
    end if
	particle(iCluster)%radiusPp = ((3d0*primParticleMatVol)/(4d0*PI))**(1d0/3d0) ! um
	particle(iCluster)%radius = particle(iCluster)%radiusPp * particle(iCluster)%nPpxP**(1d0/particle(iCluster)%fracDim) ! um

	! Sanity checks after
	if (.not. ieee_is_finite(particle(iCluster)%radius)) then
		write(*,*) 'SANITY FAIL: rp is NaN when particle is created'
		write(*,*) '  rpp, rp:', particle(iCluster)%radiusPp, particle(iCluster)%radius
		write(*,*) '  npp:', particle(iCluster)%nPpxP
		write(*,*) '  solidVolume:', particle(iCluster)%solidVolume
		write(*,*) '  mol orgC & minerals:', particle(iCluster)%molesOrgC, particle(iCluster)%molesMineral
		write(*,*) '  phase', particle(iCluster)%phase
		write(*,*) '  living', particle(iCluster)%living
		write(*,*) '  faecal', particle(iCluster)%faecal
		write(*,*) '  depth', particle(iCluster)%depth
		call WriteStatusAndStop( )
	end if
	if (particle(iCluster)%radiusPp > particle(iCluster)%radius) then
		write(*,*) 'SANITY FAIL: rpp > rp when particle is created'
		write(*,*) '  rpp, rp:', particle(iCluster)%radiusPp, particle(iCluster)%radius
		write(*,*) '  npp:', particle(iCluster)%nPpxP
		write(*,*) '  solidVolume:', particle(iCluster)%solidVolume
		write(*,*) '  mol orgC & minerals:', particle(iCluster)%molesOrgC, particle(iCluster)%molesMineral
		write(*,*) '  phase', particle(iCluster)%phase
		write(*,*) '  living', particle(iCluster)%living
		write(*,*) '  faecal', particle(iCluster)%faecal
		write(*,*) '  depth', particle(iCluster)%depth
		call WriteStatusAndStop( )
	end if
		
end subroutine ParticleRadius

! ========================================================================================

subroutine ParticlePorosity(particle, nClusters, iCluster)

	!-------------------------------------------------------------------------------------
	! Porosity values are between 0 (all solid) and 1 (all fluid). The porosity is the 
	! fraction of the aggregate volume occupied by insterstitial water, otherwise said, 
	! the fraction of an aggregate not occupied by solid matter. Only particles made by 
	! more than one primary particle have porosity.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	! Sanity checks on geometry
	if (.not. ieee_is_finite(particle(iCluster)%radiusPp) .or. particle(iCluster)%radiusPp <= 0d0) then
		write(*,*) 'ERROR: invalid radiusPp in ParticlePorosity'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(particle(iCluster)%radius) .or. particle(iCluster)%radius <= 0d0) then
		write(*,*) 'ERROR: invalid radius in ParticlePorosity'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(particle(iCluster)%fracDim)) then
		write(*,*) 'ERROR: invalid fractal dimension in ParticlePorosity'
		call WriteStatusAndStop()
	end if
	
	! Porosity calculation
	particle(iCluster)%porosity = 1d0 - (particle(iCluster)%radius &
		/particle(iCluster)%radiusPp)**(particle(iCluster)%fracDim - 3d0)

	! Sanity checks after
	if (particle(iCluster)%porosity == 1d0) then
		write(*,*) 'SANITY FAIL: P = 1 -> this particle cannot exist!'
		write(*,*) '  rpp, rp:', particle(iCluster)%radiusPp, particle(iCluster)%radius
		write(*,*) '  npp:', particle(iCluster)%nPpxP
		write(*,*) '  solidVolume:', particle(iCluster)%solidVolume
		write(*,*) '  mol orgC, TEC & minerals:', particle(iCluster)%molesOrgC, &
			particle(iCluster)%molesTepC, particle(iCluster)%molesMineral
    	call WriteStatusAndStop( )
	end if
		
end subroutine ParticlePorosity

! ========================================================================================

subroutine ParticleDensity(particle, nClusters, iCluster, waterRho)

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: waterRho
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	real*8 :: particleSolidDensity
	
	! Sanity checks
	if (.not. ieee_is_finite(particle(iCluster)%solidVolume) .or. particle(iCluster)%solidVolume <= 0d0) then
		write(*,*) 'ERROR: invalid solidVolume in ParticleDensity'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(particle(iCluster)%mass) .or. particle(iCluster)%mass <= 0d0) then
		write(*,*) 'ERROR: invalid mass in ParticleDensity'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(waterRho) .or. waterRho <= 0d0) then
		write(*,*) 'ERROR: invalid waterRho in ParticleDensity'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(particle(iCluster)%porosity) .or. &
    	particle(iCluster)%porosity < 0d0 .or. particle(iCluster)%porosity > 1d0) then
    	write(*,*) 'ERROR: invalid porosity in ParticleDensity:', particle(iCluster)%porosity
    	call WriteStatusAndStop()
	end if
	
	! Density solid material
	particleSolidDensity = (particle(iCluster)%mass/particle(iCluster)%solidVolume) * 1d12 ! g cm-3

	! Density of the aggregate (includes interstitial water)
	if (particle(iCluster)%nPpxP > 1d0) then
		particle(iCluster)%density = (1d0-particle(iCluster)%porosity)*particleSolidDensity &
			+ particle(iCluster)%porosity*waterRho ! g cm-3
	else
		particle(iCluster)%density = particleSolidDensity
	end if
	
	! Excess density necessary to calculate Stokes' sinking velocity
	particle(iCluster)%excessDensity = particle(iCluster)%density - waterRho

	! Final sanity check
	if (particle(iCluster)%density > 4d0 .or. particle(iCluster)%density < 0d0) then
		write(*,*) 'SANITY FAIL: Density calculation is wrong:', particle(iCluster)%density
		write(*,*) '  Particle type:', particle(iCluster)%initType
		write(*,*) '  Is it faecal?:', particle(iCluster)%faecal
    	call WriteStatusAndStop( )
	end if
	
end subroutine ParticleDensity

! ========================================================================================

subroutine ParticleStickiness(particle, nClusters, iCluster) 

	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	real*8 :: volTep, volSolid, volMineral, fracVolTep, fracVolMineral, S

	volSolid = particle(iCluster)%solidVolume ! um3
	volTep = 1d12*(particle(iCluster)%massTep/RHO_TEP) ! um3
	
	! Sanity check first
	if (volSolid <= 0d0 .or. .not. ieee_is_finite(volSolid)) then
		write(*,*) 'SANITY FAIL: detected wrong solid volume in stickiness calculation.'
		write(*,*) '  Mol orgC & TEC:', particle(iCluster)%molesOrgC, particle(iCluster)%molesTepC
		write(*,*) '  Mol minerals:', particle(iCluster)%molesMineral
		call WriteStatusAndStop( )
	end if
	if (volTep < 0d0 .or. .not. ieee_is_finite(volTep)) then
		write(*,*) 'SANITY FAIL: detected wrong TEP volume in stickiness calculation.'
		write(*,*) '  Mol orgC & TEC:', particle(iCluster)%molesOrgC, particle(iCluster)%molesTepC
		write(*,*) '  Mol minerals:', particle(iCluster)%molesMineral
		call WriteStatusAndStop( )
	end if
	
	! Compute volumetric fraction of TEP
  	fracVolTep = SafeDivide(volTep, volSolid, 0d0)

	! Calculate stickiness
	select case (choiceStickinessFunc)
			
	case (1) ! SLAMS-1.0 (linear scheme without mineral contribution)
		S =  MIN(1d0, MAX(0d0, fracVolTep))
	
	case (2) ! SLAMS-2.0
		! Stickiness is a saturating (Hill / Michaelis–Menten) function of TEP volume fraction, 
		! meaning low TEP additions may have little effect until a threshold, or effects saturate.
		S = SafeDivide(fracVolTep, fracVolTep + k_TEP, 0d0)
		
		! Add mineral contribution: minerals reduce adhesive contact
		volMineral = 1d12 * (particle(iCluster)%massMineral(iOpal)/RHO_OPAL &				
						   + particle(iCluster)%massMineral(iCalcite)/RHO_CALCITE &
						   + particle(iCluster)%massMineral(iClay)/RHO_CLAY) ! um3
		fracVolMineral = SafeDivide(volMineral, volSolid, 0d0)
		fracVolMineral = MIN(1d0, MAX(0d0, fracVolMineral))
		
		S = S * (1d0 - fracVolMineral)
		S = MIN(1d0, MAX(0d0, S))
		
	end select
	
	if (.not. ieee_is_finite(S)) then
		write(*,*) 'SANITY FAIL: stickiness calculation went wrong.'
		write(*,*) '  fracVolTep:', fracVolTep
		call WriteStatusAndStop( )
	end if
	
	particle(iCluster)%stickiness = S

end subroutine ParticleStickiness

! ========================================================================================

subroutine ParticleSettlingVelocity(particle, nClusters, iCluster, waterRho, waterDynVisco)

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: waterRho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle

	integer :: iVeloRecalcs 
	real*8 :: rhoW, deltaRho, signVelo, diameterP, radiusP, mu, w, wNew, reynoldsFactor, &
		reynoldsNo, dragCoeff 
	integer, parameter :: MAX_ITERS = 30   ! maximum nonlinear iterations
	real*8, parameter :: RE_MIN = 1d-12    ! floor for Reynolds number to avoid underflow
	real*8, parameter :: ABS_TOL_W = 1d-9  ! absolute tol for velocity convergence (m/s)
    real*8, parameter :: REL_TOL_W = 1d-6  ! relative tol for velocity convergence
  
	! ------ Input sanity ------
    if (.not. ieee_is_finite(waterRho) .or. waterRho < 0d0) then
        error stop 'ERROR: invalid waterRho in ParticleSettlingVelocity'
    end if
    if (.not. ieee_is_finite(waterDynVisco) .or. waterDynVisco < 0d0) then
        error stop 'ERROR: invalid waterDynVisco in ParticleSettlingVelocity'
    end if
    if (.not. ieee_is_finite(particle(iCluster)%radius) .or. particle(iCluster)%radius <= 0d0) then
        error stop 'ERROR: invalid radius in ParticleSettlingVelocity'
    end if
    
    ! ------ Unit conversions ------
	rhoW      = 1d3 * waterRho                          	! kg m-3 
	deltaRho  = ABS(particle(iCluster)%excessDensity * 1d3) ! kg m-3 
	signVelo  = SIGN(1d0, particle(iCluster)%excessDensity)
	radiusP   = particle(iCluster)%radius * 1d-6            ! m
	diameterP = 2d0 * radiusP                               ! m
	mu        = waterDynVisco * 0.1d0                   	! g cm-1 s-1 -> g m-1 s-1
	
	if (mu <= 0d0 .or. deltaRho == 0d0) then
	  	particle(iCluster)%velocity = 0d0
	  	return
	end if

	! ------ Stokes initial guess ------
	w = (2d0/9d0) * GRAVITY_CNT * deltaRho * (radiusP*radiusP) / mu ! m s-1
	particle(iCluster)%velocity = w * SECONDS_PER_DAY ! m d-1 (update to insert into Reynolds number calculation)

	! ------ Calculate Reynolds number ------
	call CalculateParticleReynoldsNumber(reynoldsFactor, reynoldsNo, particle, nClusters, iCluster, waterRho, waterDynVisco)
	
	! ------ Modify initial guess according to Re regime ------
	
	! If near neutral flow (i.e. small Re) --> stick with Stokes 
	if (reynoldsNo < REYNOLDS_LAMINAR_LIMIT) then
	  	particle(iCluster)%velocity = signVelo * w * SECONDS_PER_DAY ! m d-1
	  	return
	end if
	
	! Otherwise iterate with White's Cd(Re), and update w as you update reynoldsNo. Stop when 
	! w converges within combined absolute/relative tol
	do iVeloRecalcs = 1, MAX_ITERS
		reynoldsNo = MAX(reynoldsNo, RE_MIN)
		
		! Drag coefficient
	  	dragCoeff  = (24d0/reynoldsNo) + (6d0/(1d0 + SQRT(reynoldsNo))) + 0.4d0
	  	if (.not. ieee_is_finite(dragCoeff) .or. dragCoeff <= 0d0) exit
	  	
	  	! Velocity formula
	  	wNew = (4d0 * deltaRho * diameterP * GRAVITY_CNT) / (3d0 * rhoW * dragCoeff)
        if (wNew <= 0d0) exit
        wNew = SQRT(wNew) ! m s-1
        if (.not. ieee_is_finite(wNew)) exit
        if (ABS(wNew - w) <= MAX(ABS_TOL_W, REL_TOL_W * wNew)) exit ! tolerance check
	  	w = wNew
	  	
	  	! Re
	  	reynoldsNo = MAX(ABS(w) * reynoldsFactor, RE_MIN) ! reynoldsFactor does not change with w as it does not depend on 'w'
	end do
	
	! ------ Restore physical sign: + sinking, − rising ------
	particle(iCluster)%velocity = signVelo*w*SECONDS_PER_DAY ! m d-1
	
	! ------ Final sanity check ------
	if (.not. ieee_is_finite(particle(iCluster)%velocity)) then
		call PrintParticleProperties(particle, nClusters, iCluster, 'ERROR: ParticleSettlingVelocity gone wrong')
    	call WriteStatusAndStop( )
	end if
	
end subroutine ParticleSettlingVelocity

! ========================================================================================

subroutine CalculateParticleReynoldsNumber(ReFactor, ReNumber, particle, nClusters, iCluster, &
	waterRho, waterDynVisco)
	
	!-------------------------------------------------------------------------------------
	! Calculate the hydrodynamic Reynolds number of the particle based on the particle's
	! sinking velocity, diameter, and water viscosity. It tells you whether the flow around 
	! the particle itself is laminar or inertial because of the particle's motion through 
	! the water.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: waterRho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	real*8, intent(out) :: ReFactor, ReNumber 
	
	real*8 :: rhoW, diameterP, radiusP, mu, w
	real*8, parameter :: RE_MIN = 1d-12 ! floor for Reynolds number to avoid underflow
	
	! ------ Initialise these (default output) ------
	ReFactor = 0d0
  	ReNumber = 0d0
  	
	! ------ Input sanity ------
    if (.not. ieee_is_finite(particle(iCluster)%radius) .or. particle(iCluster)%radius <= 0d0) then
    	write(*,*) 'ERROR: invalid radius in CalculateParticleReynoldsNumber'
		call WriteStatusAndStop()
    end if
    if (.not. ieee_is_finite(particle(iCluster)%velocity)) then
    	write(*,*) 'ERROR: invalid velocity in CalculateParticleReynoldsNumber'
		call WriteStatusAndStop()
    end if

	! ------ Unit conversion ------
	rhoW      = waterRho * 1d3                              ! g cm-3 --> kg m-3                        
	radiusP   = particle(iCluster)%radius * 1d-6            ! um --> m
	diameterP = 2d0 * radiusP                               ! m
	mu        = waterDynVisco * 0.1d0                       ! g cm-1 s-1 -> kg m-1 s-1
	w         = particle(iCluster)%velocity/SECONDS_PER_DAY ! m d-1 --> m s-1
	
	! ------ Reynolds number calculations ------
	ReFactor = rhoW * diameterP / mu
	ReNumber = MAX(ABS(w) * ReFactor, RE_MIN) ! dimensionless
	
end subroutine CalculateParticleReynoldsNumber

! ========================================================================================

subroutine CalculateTurbulenceScaleReynoldsNumber(ReTurb, particle, nClusters, iCluster, &
	kolmogorovLengthScale)
	
    !-------------------------------------------------------------------------------------
    ! Calculate turbulence-scale Reynolds number associated with eddies of size comparable 
    ! to the particle. It comes from Kolmogorov / viscous-subrange scaling and measures 
    ! the intensity of turbulent eddies at the particle scale (i.e., whether turbulent 
    ! eddies acting on the particle are energetic enough to break it up).
    !-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: kolmogorovLengthScale
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	real*8, intent(out) :: ReTurb 
	
	real*8 :: radiusP
	
	! ------ Initialise these (default output) ------
	ReTurb = 0d0
	
	! ------ Input sanity ------
    if (.not. ieee_is_finite(particle(iCluster)%radius) .or. particle(iCluster)%radius <= 0d0) then
    	write(*,*) 'ERROR: invalid radius in CalculateTurbulenceScaleReynoldsNumber'
		call WriteStatusAndStop()
    end if
    if (kolmogorovLengthScale == 0d0) return ! no turbulence if Kolmogorov scale is zero
	
	! ------ Unit conversion ------
	radiusP = particle(iCluster)%radius * 1d-6 ! um --> m
	
	! ---------- Turbulence-scale Reynolds number ----------
	if (radiusP <= kolmogorovLengthScale) then
		! Viscous subrange: use ReTurb ≈ 4 * (r/eta)^2
		ReTurb = 4d0 * (radiusP/kolmogorovLengthScale)**2
	else
		! Inertial range: adopt approximate scaling ∝ (r/eta)^(4/3)
		! choose prefactor consistent with viscous value at r=eta (4.0)
		ReTurb = 4d0 * (radiusP/kolmogorovLengthScale)**(4d0/3d0)
	end if
	
	! ---------- Final sanity ----------
    if (.not. ieee_is_finite(ReTurb) .or. ReTurb < 0d0) then
        write(*,*) 'ERROR: ReTurb invalid in CalculateTurbulenceScaleReynoldsNumber'
        call WriteStatusAndStop()
    end if

end subroutine CalculateTurbulenceScaleReynoldsNumber

! ========================================================================================

subroutine ComputeParticleAttributes(particle, nClusters, iCluster, particleInitType, &
	iPft, isLiving, isFaecal, iTimeStep, OmegaCalc, waterDepth, nParticlesPerCluster, &
	carbonQuota, waterRho, waterDynVisco)

	!-------------------------------------------------------------------------------------
	! This subroutine is used only at the time of particle creation for phytoplankton cells, 
	! TEP particles, clay particles, and mesozooplankton dead bodies. It is not used for 
	! faecal pellets, aggregates, or for any future time steps after the initial creation 
	! of the particle.
	!-------------------------------------------------------------------------------------
		
	integer, intent(in) :: nClusters, iCluster, particleInitType, iPft, isLiving, &
		isFaecal, iTimeStep
	real*8, intent(in) :: OmegaCalc, waterDepth, nParticlesPerCluster, carbonQuota, &
		waterRho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle

	real*8 :: adjustedOmegaCalc
	
	! ------ Sanity checks first ------
	if (nParticlesPerCluster < 1d0) then
		write(*,*) 'ERROR: nParticlesPerCluster < 1 when particle is created'
		call WriteStatusAndStop( )
	end if
	if (.not. ieee_is_finite(waterDepth) .or. waterDepth < 0d0) then
		write(*,*) 'ERROR: particle depth is NaN when particle is created'
		call WriteStatusAndStop( )
	end if
	if (.not. ieee_is_finite(carbonQuota) .or. carbonQuota < 0d0) then
		write(*,*) 'ERROR: invalid carbonQuota in ComputeParticleAttributes'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(waterRho) .or. waterRho <= 0d0) then
		write(*,*) 'ERROR: invalid waterRho in ComputeParticleAttributes'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(waterDynVisco) .or. waterDynVisco <= 0d0) then
		write(*,*) 'ERROR: invalid waterDynVisco in ComputeParticleAttributes'
		call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(OmegaCalc)) then
		write(*,*) 'ERROR: invalid OmegaCalc in ComputeParticleAttributes'
		call WriteStatusAndStop()
	end if
	if (particleInitType < 1 .or. particleInitType > 7) then
		write(*,*) 'ERROR: invalid particleInitType in ComputeParticleAttributes:', particleInitType
		call WriteStatusAndStop()
	end if
	if (isLiving /= 0 .and. isLiving /= 1) then
		write(*,*) 'ERROR: isLiving must be 0 or 1 in ComputeParticleAttributes'
		call WriteStatusAndStop()
	end if
	if (isFaecal /= 0 .and. isFaecal /= 1) then
		write(*,*) 'ERROR: isFaecal must be 0 or 1 in ComputeParticleAttributes'
		call WriteStatusAndStop()
	end if

	! ------ Particle settings ------
	particle(iCluster)%id = iCluster
	particle(iCluster)%initType = particleInitType
	particle(iCluster)%phase = 1	
	particle(iCluster)%living = isLiving
	particle(iCluster)%faecal = isFaecal
	particle(iCluster)%tstepCreat = iTimeStep		
	particle(iCluster)%depth = waterDepth ! m

	! ------ Particle counts ------
	particle(iCluster)%nPpxC = nParticlesPerCluster
	particle(iCluster)%nPpxP = 1d0 ! means this is not an aggregate particle
	particle(iCluster)%nPxC = nParticlesPerCluster
	
	! ------ Initialise material content ------
	particle(iCluster)%molesOrgC = 0d0
	particle(iCluster)%massOrgMatter = 0d0
	particle(iCluster)%molesTepC = 0d0
	particle(iCluster)%massTep = 0d0
	particle(iCluster)%molesMineral(:) = 0d0
	particle(iCluster)%massMineral(:) = 0d0

	! ------ Compute material content according to type ------
	! Phytoplankton cell
	if (isLiving == 1) then 
		if (iPft < 1 .or. iPft > nPfts) then
			write(*,*) 'ERROR: invalid iPft in ComputeParticleAttributes:', iPft
			call WriteStatusAndStop()
		end if
		particle(iCluster)%molesOrgC = carbonQuota
		particle(iCluster)%stickiness = 0d0 ! 0-1
		particle(iCluster)%initPft = iPft
		
		if (iPft == 1) then 		! diatom: silicification
			particle(iCluster)%molesMineral(iOpal) = carbonQuota*Si2C_diat
			
		else if (iPft == 3) then 	! cocco: calcification
    		adjustedOmegaCalc = MAX(OmegaCalc, 1.01d0) ! ensure OmegaCalc is above 1 for calcification (slightly above 1 to avoid division by zero)
			particle(iCluster)%molesMineral(iCalcite) = carbonQuota * Calc2C_cocco_max &
				* ((adjustedOmegaCalc-1d0)/(k_omega + (adjustedOmegaCalc-1d0)))
			
		end if
		
		! Notice we give phytoplankton cells an initial stickiness
		particle(iCluster)%stickiness = stickiness_phyto_cell_initial ! 0-1

	! TEP											
	else if (particleInitType == 5) then	
		particle(iCluster)%molesTepC = carbonQuota
		particle(iCluster)%stickiness = 1d0 ! 0-1

	! Terrigenous particle (clay)						
	else if (particleInitType == 6) then 			
		particle(iCluster)%molesMineral(iClay) = clay_quota
		particle(iCluster)%stickiness = 0d0 ! 0-1

	! Zooplankton dead body						
	else if (particleInitType == 7) then
		particle(iCluster)%molesOrgC = carbonQuota	
		particle(iCluster)%stickiness = 0d0 ! 0-1	
	
	end if

	! ------ Compute attributes ------
	call ParticleFractalDimension(particle, nClusters, iCluster, particleInitType)
	call ParticleDryMass(particle, nClusters, iCluster)
	call ParticleMaterialVolume(particle, nClusters, iCluster) 
	call ParticleRadius(particle, nClusters, iCluster) 
	call ParticlePorosity(particle, nClusters, iCluster)
	call ParticleDensity(particle, nClusters, iCluster, waterRho) 
	call ParticleSettlingVelocity(particle, nClusters, iCluster, waterRho, waterDynVisco) 

	! ------ Sanity check ------
	call CheckParticleSanity(particle, nClusters, iCluster, 'primary particle creation')
		
end subroutine ComputeParticleAttributes

! ========================================================================================

subroutine CompactParticleArray(particle, nClusters, iLastLocus, nPhase1clusters, phase1clusterIndices)

	!-------------------------------------------------------------------------------------
	! This subroutine squeezes the particle array by getting rid of clusters that are not 
	! in phase 1. It is used after first identifying those clusters.
	!-------------------------------------------------------------------------------------
	
	integer, intent(in) :: nClusters, nPhase1clusters
	integer, dimension(nPhase1clusters), intent(in) :: phase1clusterIndices
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	integer, intent(inout) :: iLastLocus
	
	integer :: iSrc, iDst, iClean
	type(lagrangianStateVars) :: tmpParticle

	! ------ Compact phase-1 clusters ------
	do iDst = 1, nPhase1clusters ! destination index
		iSrc = phase1clusterIndices(iDst) ! source index
		if (iSrc < 1 .or. iSrc > iLastLocus) then
        	write(*,*) 'ERROR: invalid phase1clusterIndices entry'
        	write(*,*) '  index=', iDst, ' value=', iSrc
        	write(*,*) '  iLastLocus=', iLastLocus
        	call WriteStatusAndStop()
    	end if
		tmpParticle = particle(iSrc)
    	particle(iDst) = tmpParticle ! write back to compacted slot
    	particle(iDst)%id = iDst ! ensure id matches new position
	end do
	
	! ------ Clean trailing slots ------
	do iClean = nPhase1clusters+1, iLastLocus
        call InitialiseParticle(particle(iClean))
	end do

	! ------ Update iLastLocus to reflect new end of active region ------
	iLastLocus = nPhase1clusters

end subroutine CompactParticleArray

! ========================================================================================

subroutine ShrinkParticle(particle, nClusters, iCluster, volMaterialDegraded)

	!-------------------------------------------------------------------------------------
	! This subroutine applies to aggregates that initially had more than two primary particles,
	! but that have undergone respiration, solubilisation or mineral dissolution. After material 
	! is transferred to the dissolved pool in these processes, the remaining particulate material 
	! is redistributed in a smaller number of primary particles. These processes can 
	! end up proceeding far enough that an aggregate becomes a single primary particle 
	! during the shrinking process.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: volMaterialDegraded
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	integer :: nDegradedPrimParticles
	real*8 :: nCurrPrimParticles, volPrimParticle, expectedLoss, maxLoss, baseLoss, frac, &
		harvest, newPrimaryParticleCount

	! ------ Sanity checks first ------
  	if (.not. ieee_is_finite(volMaterialDegraded) .or. volMaterialDegraded < 0d0) then
   		write(*,*) 'ERROR: invalid volMaterialDegraded in ShrinkParticle'
   		call WriteStatusAndStop( )
	end if
	if (.not. ieee_is_finite(particle(iCluster)%solidVolume) .or. particle(iCluster)%solidVolume <= 0d0) then
   		write(*,*) 'ERROR: invalid solid volume in ShrinkParticle'
   		call WriteStatusAndStop( )
	end if
	if (.not. ieee_is_finite(particle(iCluster)%nPpxP) .or. particle(iCluster)%nPpxP < 1d0) then
   		write(*,*) 'ERROR: invalid nPpxP in ShrinkParticle'
   		call WriteStatusAndStop( )
	end if

	nCurrPrimParticles = particle(iCluster)%nPpxP
	volPrimParticle = particle(iCluster)%solidVolume / nCurrPrimParticles ! um3
	if (.not. ieee_is_finite(volPrimParticle) .or. volPrimParticle <= 0d0) then
	   	write(*,*) 'ERROR: invalid volPrimParticle in ShrinkParticle'
   		call WriteStatusAndStop( )
	end if
	
	! ---------- Expected loss ----------
  	expectedLoss = volMaterialDegraded / volPrimParticle ! how many whole monomers does that degraded volume represent?
  	if (.not. ieee_is_finite(expectedLoss) .or. expectedLoss <= 0d0) return
  	
  	! ---------- Maximum allowed loss ----------
  	maxLoss = nCurrPrimParticles - 1d0 ! max. monomers you’re allowed to remove so one monomer always remains
  	if (maxLoss <= 0d0) return
  	
  	expectedLoss = MIN(expectedLoss, maxLoss)
  	
  	! ------ Stochastic rounding ------
  	baseLoss = SafeFloorNonNegative(expectedLoss)
    frac = expectedLoss - baseLoss
    frac = MAX(0d0, MIN(1d0, frac))
    call random_number(harvest)
    if (harvest < frac) then
        nDegradedPrimParticles = baseLoss + 1d0
    else
        nDegradedPrimParticles = baseLoss
    end if
    
    ! ------ Safety clamp -------
    if (nDegradedPrimParticles < 0d0) nDegradedPrimParticles = 0d0
    if (nDegradedPrimParticles > maxLoss) nDegradedPrimParticles = maxLoss

	! ------ Update particle counts ------
	newPrimaryParticleCount = nCurrPrimParticles - nDegradedPrimParticles
	if (.not. ieee_is_finite(newPrimaryParticleCount) .or. newPrimaryParticleCount < 1d0) then
        write(*,*) 'ERROR: invalid newPrimaryParticleCount in ShrinkParticle:', newPrimaryParticleCount
    	call WriteStatusAndStop()
    end if
	particle(iCluster)%nPpxP = newPrimaryParticleCount
	particle(iCluster)%nPpxC = particle(iCluster)%nPpxP * particle(iCluster)%nPxC
	
	! ------ Update particle type ------
	particle(iCluster)%initType = ParticleUpdatedType(particle, nClusters, iCluster)
		
end subroutine ShrinkParticle			

! ========================================================================================

function VolumetricFractionOfOrganicMatter(particle, nClusters, iCluster) result(vfrac)

	!-------------------------------------------------------------------------------------
	! VolumetricFractionOfOrganicMatter returns the volume fraction of non-TEP organic 
	! matter when present, and uses TEP matter only as a proxy for organic material when 
	! non-TEP OM is absent
	!-------------------------------------------------------------------------------------
	
	real*8 :: vfrac
	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	real*8 :: volMineral, volOrganicMatter, den
	
	! ------ Initialise these (default output) ------
	vfrac = 0d0
	
	! ------ Input sanity ------
	if (.not. ieee_is_finite(particle(iCluster)%massMineral(iOpal)) .or. &
        .not. ieee_is_finite(particle(iCluster)%massMineral(iCalcite)) .or. &
        .not. ieee_is_finite(particle(iCluster)%massMineral(iClay))) then
        write(*,*) 'ERROR: invalid mass mineral in VolumetricFractionOfOrganicMatter'
        call WriteStatusAndStop()
    end if
    if (.not. ieee_is_finite(particle(iCluster)%massOrgMatter) .or. &
        .not. ieee_is_finite(particle(iCluster)%massTep)) then
        write(*,*) 'ERROR: invalid mass organic in VolumetricFractionOfOrganicMatter'
        call WriteStatusAndStop()
    end if
    
    ! ---------- Mineral volume ----------
    volMineral = 0d0 ! cm3
    if (particle(iCluster)%massMineral(iOpal)    > 0d0) volMineral = volMineral + particle(iCluster)%massMineral(iOpal)/RHO_OPAL
    if (particle(iCluster)%massMineral(iCalcite) > 0d0) volMineral = volMineral + particle(iCluster)%massMineral(iCalcite)/RHO_CALCITE
    if (particle(iCluster)%massMineral(iClay)    > 0d0) volMineral = volMineral + particle(iCluster)%massMineral(iClay)/RHO_CLAY

	! ---------- Organic volume ----------
	volOrganicMatter = 0d0
	
	! Prefer OM (excluding TEP); fallback to TEP only if OM is zero	
	if (particle(iCluster)%massOrgMatter == 0d0 .and. particle(iCluster)%massTep > 0d0) then
		volOrganicMatter = particle(iCluster)%massTep/RHO_TEP ! cm3
	else if (particle(iCluster)%massOrgMatter > 0d0) then ! exclude TEP
		volOrganicMatter = particle(iCluster)%massOrgMatter/RHO_ORGMATTER ! cm3
	end if
	
	! ------ Fraction ------
	den = volOrganicMatter + volMineral
	if (.not. ieee_is_finite(den) .or. den <= 0d0) return
  	vfrac = volOrganicMatter / den
    vfrac = MIN(1d0, MAX(0d0, vfrac)) ! clamp [0, 1]
    				
end function VolumetricFractionOfOrganicMatter

! ========================================================================================

function FaecalPelletPorosity(particle, nClusters, iCluster) result(fppor)

	!-------------------------------------------------------------------------------------
	! The porosity of a faecal pellet is calculated differently from any other aggregate.
	!-------------------------------------------------------------------------------------
	
	real*8 :: fppor
	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	real*8 :: fracVolOrgMatter, fracVolMineral
	
	! ------ Initialise these (default output) ------
	fppor = porosity_faecalpell_min
	
	! ------ Calculate fractions ------
	fracVolOrgMatter = VolumetricFractionOfOrganicMatter(particle, nClusters, iCluster)
	if (.not. ieee_is_finite(fracVolOrgMatter)) then
        write(*,*) 'ERROR: invalid fracVolOrgMatter in FaecalPelletPorosity'
        call WriteStatusAndStop()
    end if
    fracVolOrgMatter = MIN(1d0, MAX(0d0, fracVolOrgMatter)) ! ensure its physically meaningful
	fracVolMineral = 1d0 - fracVolOrgMatter 
	
	! ------ Porosity calculation ------
	select case (choiceMineralEffOnFaecalPellPorosity)
	case (1) ! minerals increase porosity
		fppor = porosity_faecalpell_min + &
			(porosity_faecalpell_max - porosity_faecalpell_min)*fracVolMineral
	case (2) ! minerals decrease porosity
		fppor = porosity_faecalpell_min + &
			(porosity_faecalpell_max - porosity_faecalpell_min)*(1d0-fracVolMineral)
	end select
	fppor = MIN(1d0, MAX(0d0, fppor)) ! clamp [0, 1]
	
end function FaecalPelletPorosity

! ========================================================================================

function ParticleUpdatedType(particle, nClusters, iCluster) result(updatedtype)

	!-------------------------------------------------------------------------------------
    ! Update particle type for non-living particles after processes that may alter
    ! particle's composition or structure, namely:
    !   - coagulation and breakup,
    !   - particle shrinking due to respiration, solubilisation, or mineral dissolution.
    !
    ! The function assumes that living phytoplankton cells have already been excluded and
    ! should therefore never be passed to this function.
    !-------------------------------------------------------------------------------------

	integer :: updatedtype
	integer, intent(in) :: nClusters, iCluster
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	integer :: isFaecal, isLiving, particleInitType
	real*8 :: nPrimParticlesPerParticle, totalOrgC, totalMineral
	logical :: isAggregate, isSingleLiving, isSingleDead

	updatedtype = 0
	
	nPrimParticlesPerParticle = particle(iCluster)%nPpxP 
	isFaecal = particle(iCluster)%faecal
	isLiving = particle(iCluster)%living
	particleInitType = particle(iCluster)%initType
	
    if (isLiving == 1) then
        write(*,*) 'ERROR: trying to update particle type of a living phytoplankton cell in ParticleUpdatedType'
        call WriteStatusAndStop()
    end if
	
!-----------------------------------------------------------------------------------------
! These checks enforce logical consistency of particle attributes but significantly increase 
! runtime when called repeatedly within SLAMS. They are therefore omitted in production 
! runs, as the corresponding properties are validated upstream.
!-----------------------------------------------------------------------------------------
!	
! 	! Validate inputs
!     if (isFaecal < 0 .or. isFaecal > 1) then
!         print *, "ERROR: isFaecal must be 0 or 1."
!         open(unit=10, file="program_status.tmp", status="replace", action="write")
!     	write(10, '(A)') "error"
!     	close(10)
! 		stop
!     end if
!     if (isLiving < 0 .or. isLiving > 1) then
!         print *, "ERROR: isLiving must be 0 or 1."
!         open(unit=10, file="program_status.tmp", status="replace", action="write")
!     	write(10, '(A)') "error"
!     	close(10)
! 		stop
!     end if
!     if (particleInitType < 1 .or. particleInitType > 9) then
!         print *, "ERROR: particleInitType must be between 1 and 9."
!         open(unit=10, file="program_status.tmp", status="replace", action="write")
!     	write(10, '(A)') "error"
!     	close(10)
! 		stop
!     end if
!     if (.not. ieee_is_finite(nPrimParticlesPerParticle) .or. nPrimParticlesPerParticle < 1d0) then
!         write(*,*) 'ERROR: invalid nPrimParticlesPerParticle in ParticleUpdatedType:', nPrimParticlesPerParticle
!         open(unit=10, file="program_status.tmp", status="replace", action="write")
!     	write(10, '(A)') "error"
!     	close(10)
! 		stop
!     end if

    totalOrgC = particle(iCluster)%molesOrgC + particle(iCluster)%molesTepC
    totalMineral = particle(iCluster)%massMineral(iOpal) + particle(iCluster)%massMineral(iCalcite)
	if (isFaecal == 1) then
		updatedtype = 1	 ! faecal particle
		return
	elseif (isFaecal == 0 .and. nPrimParticlesPerParticle > 1d0) then
		updatedtype = 2  ! non-faecal aggregate
		return
	elseif (isFaecal == 0 .and. nPrimParticlesPerParticle == 1d0 .and. totalOrgC > 0d0) then
		updatedtype = 9 ! organic single primary particle fragment
		return
	elseif (isFaecal == 0 .and. nPrimParticlesPerParticle == 1d0 .and. totalOrgC == 0d0 &
		.and. totalMineral > 0d0 .and. particle(iCluster)%massMineral(iClay) == 0d0) then
		updatedtype = 8 ! inorganic single primary particle fragment
		return
	elseif (isFaecal == 0 .and. nPrimParticlesPerParticle == 1d0 .and. totalOrgC == 0d0 &
		.and. totalMineral == 0d0 .and. particle(iCluster)%massMineral(iClay) > 0d0) then
		updatedtype = 6  ! clay fragment
		return	
	else
		write(*,*) 'ERROR: invalid ParticleUpdatedType combination:', updatedtype
		call WriteStatusAndStop()	
	end if
	
end function ParticleUpdatedType

! ========================================================================================

end module calcparticleattributes
