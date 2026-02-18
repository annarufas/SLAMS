#include "blockdefinitions.h"

module particledynamics

! ----------------------------------------------------------------------------------------
! This module handles the collision and aggregation of particles, fragmenting them, 
! packing clusters and disaggregating clusters.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use safemath, only: SafeFloorNonNegative, SafeDivide, SafeLog, SafeCeilingNonNegative, &
	SafeSplitIntegerLike
use modelcounters, only: nTimesEnteringCollLoop, nCollisionsType1, nCollisionsType2, &
	nCollisionsType3, nCollisionsType4, nFragmentedBigClusters 
use modelconstants, only: SECONDS_PER_DAY, BOLTZMANN_CNT, PI
use modelparameters, only: maxNumAuxTerms, maxNumSmsTerms, timeStep, choiceCollisionKernels, &
	choiceIncludeHydrodynamicForces, iCalcite, iOpal, iClay, nMinerals, iMicrobSolubOrgC, &
	iCoagulationSuccess, iCollisionKernel, iNumClusterPairsEvalCoagu, iBrownianKernel, &
	iShearKernel, iSettlingKernel, iCoagulationProbability, iNumPxCmore, iNumPxCless, &
	operational_size_poc_min, detection_limit_poc, choiceCollisionSamplingScheme, &
	breaking_reynolds_threshold, choiceCriteriaToBreak, decreaseCollisionTimeframeFactor, &
	maxPrimParticlesPerAggregateForCollision
use sanitychecks, only: CheckParticleSanity, WriteStatusAndStop, NearlyEqual, IsExceedingInitialAmount
use calcparticleattributes, only: ParticleFractalDimension, ParticleDryMass, ParticleMaterialVolume, &
	ParticleRadius, ParticlePorosity, ParticleDensity, ParticleStickiness, ParticleSettlingVelocity, &
	VolumetricFractionOfOrganicMatter, ParticleUpdatedType, CalculateParticleReynoldsNumber, &
	CalculateTurbulenceScaleReynoldsNumber
		
implicit none
private
public :: CoagulateParticles, BreakUnstableParticles, BreakParticle

contains

! ========================================================================================

subroutine CoagulateParticles(particle, nClusters, nLayerClusters, layerClusterIndices, &
	SMSterm, auxTerm, auxCount, TempK, Rho, waterDynVisco, waterShearRate, kolmogorovLengthScale, &
	gridCellVolume, iProfile)

	integer, intent(in) :: nClusters, nLayerClusters, iProfile
	integer, dimension(nLayerClusters), intent(in) :: layerClusterIndices	
	real*8, intent(in) :: TempK, Rho, waterDynVisco, waterShearRate, kolmogorovLengthScale, gridCellVolume
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle	
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm	
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount

	integer :: iCollision, iClusterA, iClusterB, nClustersElig, nTotalPairs, nPairsSubset, &
		posA, posB, tmp
	real*8 :: harvest
	integer, allocatable :: eligParticles(:)
	integer, allocatable :: perm(:)
	integer :: i, idx, jidx

	! Build an eligible list of particles that will be evaluated for collision
	if (nLayerClusters < 2) return
	allocate(eligParticles(nLayerClusters))
	nClustersElig = 0
	do i = 1, nLayerClusters
	  	idx = layerClusterIndices(i)
	  	if (IsClusterIndividuallyEligibleForCollision(particle, nClusters, idx, Rho, waterDynVisco, kolmogorovLengthScale)) then
			nClustersElig = nClustersElig + 1
			eligParticles(nClustersElig) = idx
 	  	end if
	end do
	if (nClustersElig < 2) then
	  	deallocate(eligParticles)
	  	return
	end if
	
	! Count number of unordered pairs of clusters (i.e., combinations). In combinations, 
	! we are NOT allowed to count the same pair again, regardless of the order 
	! (e.g., A and B = B and A), that's why we divide by 2
	nTotalPairs = nClustersElig*(nClustersElig-1)/2
	
	select case (choiceCollisionSamplingScheme)
	
	!-------------------------------------------------------------------------------------
	! OPTION 1: As in SLAMS-1.0
	! Stochastically samples all combinatorial pairs, and there is replacement (i.e., a 
	! cluster can be picked many times across iterations). This method is stochastic, so
	! results vary run-to-run.
	! Not correct if our aim is to evaluate each unique pair once or to perform SDM-style 
	! disjoint pairing. Also inefficient: doing nTotalPairs draws with replacement will 
	! generally not cover all unique pairs and will repeat pairs many time.
	!-------------------------------------------------------------------------------------
	
	case (1)
		nPairsSubset = nTotalPairs ! since there is no subsetting (it is expensive)
		
		do iCollision = 1, nTotalPairs	
		
			! Draw posA uniformly in 1..nClustersElig  
			call random_number(harvest); posA = 1 + INT(harvest * nClustersElig)     
		
			! Draw posB uniformly from the n-1 indices != posA:
			call random_number(harvest); posB = 1 + INT(harvest * (nClustersElig-1))
			if (posB >= posA) posB = posB + 1   ! map into 1..nClustersElig but != posA
			
			! Order them the pair explicitly as an unordered pair (so A,B and B,A do not 
			! appear in different draws)
			if (posA > posB) then
				tmp = posA; posA = posB; posB = tmp
			end if
		
			iClusterA = eligParticles(posA)
			iClusterB = eligParticles(posB)
			
			! Guards
			if (ShouldSkipCollision(particle, nClusters, iClusterA, iClusterB, Rho, waterDynVisco, kolmogorovLengthScale)) cycle
			call CollideParticles(particle, nClusters, iClusterA, iClusterB, nTotalPairs, &
				nPairsSubset, SMSterm, auxTerm, auxCount, TempK, Rho, waterDynVisco, &
				waterShearRate, kolmogorovLengthScale, gridCellVolume, iProfile)
		end do

	!-------------------------------------------------------------------------------------
	! OPTION 2: samples all particles, and every unique unordered pair {A<B} is evaluated 
	! exactly once (i.e., sampling with no replacement). This method is deterministic and 
	! exact (samples every pair). Combinatorial enumeration instead of MC sampling. Cost 
	! is O(N²), but no sampling bias.
	!-------------------------------------------------------------------------------------
	
	case (2)
		nPairsSubset = nTotalPairs ! since there is no subsetting (it is expensive)

		do posA = 1, nClustersElig-1
			do posB = posA + 1, nClustersElig
				iClusterA = eligParticles(posA)
				iClusterB = eligParticles(posB)
		
				! Guards 
				if (iClusterA == iClusterB) cycle
				if (ShouldSkipCollision(particle, nClusters, iClusterA, iClusterB, Rho, waterDynVisco, kolmogorovLengthScale)) cycle
				call CollideParticles(particle, nClusters, iClusterA, iClusterB, nTotalPairs, &
					nPairsSubset, SMSterm, auxTerm, auxCount, TempK, Rho, waterDynVisco, &
					waterShearRate, kolmogorovLengthScale, gridCellVolume, iProfile)
			end do
		end do

	!-------------------------------------------------------------------------------------
	! OPTION 3: as in the Super-Droplet Method, where there is stochastic subsampling (O(n) cost), 
	! and sampling without replacement (i.e., no droplet used twice).
	! Take nPairsSubset = nClustersElig/2 disjoint pairs 
	!-------------------------------------------------------------------------------------
	
	case (3)
	
		! Calculate subsample size so that computational cost scales as O(n), i.e., linear 
		! in the number of super-droplets
		nPairsSubset = nClustersElig / 2 
		if (nPairsSubset < 1) then
			deallocate(eligParticles); return
		end if
		
		! Allocate permutation array (1..nClustersElig)
		allocate(perm(nClustersElig))
		do i = 1, nClustersElig
			perm(i) = i
		end do
		
		! Fisher-Yates shuffle to randomise eligible indices
		do i = nClustersElig, 2, -1
			call random_number(harvest)                 ! harvest in [0,1)
			jidx = INT(harvest * dble(i)) + 1           ! uniform in 1..i
			if (jidx < 1) jidx = 1
			if (jidx > i) jidx = i
			tmp = perm(i); perm(i) = perm(jidx); perm(jidx) = tmp
		end do
	
		! Now take first 2*nPairsSubset entries and form pairs (perm(1),perm(2)), (perm(3),perm(4)), ...
		do iCollision = 1, nPairsSubset
			posA = perm(2*iCollision - 1)    ! position in 1..nClustersElig
			posB = perm(2*iCollision)        ! position in 1..nClustersElig
	
			! Map permutation positions -> actual cluster indices
			iClusterA = eligParticles(posA)
			iClusterB = eligParticles(posB)
			
			! Guards 
			if (iClusterA == iClusterB) cycle
			if (ShouldSkipCollision(particle, nClusters, iClusterA, iClusterB, Rho, waterDynVisco, kolmogorovLengthScale)) cycle
			call CollideParticles(particle, nClusters, iClusterA, iClusterB, nTotalPairs, &
				nPairsSubset, SMSterm, auxTerm, auxCount, TempK, Rho, waterDynVisco, &
				waterShearRate, kolmogorovLengthScale, gridCellVolume, iProfile)
		end do
		deallocate(perm)
		
	end select
	
	auxTerm(iNumClusterPairsEvalCoagu) = auxTerm(iNumClusterPairsEvalCoagu) + dble(nTotalPairs)
	auxCount(iNumClusterPairsEvalCoagu) = auxCount(iNumClusterPairsEvalCoagu) + 1
	deallocate(eligParticles)
	
end subroutine CoagulateParticles

! ========================================================================================

logical function IsClusterIndividuallyEligibleForCollision(particle, nClusters, iCluster, &
	Rho, waterDynVisco, kolmogorovLengthScale) result(eligible)

	!-------------------------------------------------------------------------------------
    ! Determine whether an individual particle cluster is eligible to participate in
    ! coagulation.
    !
    ! This helper function is used by CoagulateParticles and applies the same exclusion
    ! criteria as ShouldSkipCollision below, but at the level of a single cluster.
    !
    ! NOTE ON THE HARD-WIRED nPpxP UPPER BOUND (1e9):
    ! An upper limit of 1e9 primary particles per particle (nPpxP) is enforced here, as
    ! inherited from SLAMS-1.0. The origin of this threshold is currently ambiguous:
    ! it may represent either (i) a modelling contract or (ii) a numerical limitation.
    !
    ! If interpreted as a numerical limitation, increasing this cap causes the particle
    ! population to become numerically unstable and the particle landscape to deteriorate.
    ! If interpreted as a modelling contract, clusters exceeding this value are considered
    ! physically unstable and are excluded from further coagulation, as they would more
    ! plausibly fragment rather than continue to grow.
    !
    ! The value is therefore treated as a global assumption throughout SLAMS and must not
    ! be modified locally. In future developments, this hard threshold could be replaced
    ! by a physically based stability criterion (e.g. Reynolds-regime-dependent breakup;
    ! see commented code below).
    !-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: Rho, waterDynVisco, kolmogorovLengthScale
  	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle

  	eligible = .false.

  	! Basic per-cluster acceptance tests:
  	if (particle(iCluster)%phase /= 1) return
  	if (particle(iCluster)%nPxC < 1d0) return ! choose >=1, that lets single-particle clusters appear in the eligible set (useful if you allow single+multi collisions)
  	if (2d0*particle(iCluster)%radius < operational_size_poc_min) return
  	if (particle(iCluster)%nPpxP >= maxPrimParticlesPerAggregateForCollision) return

!   Optional stability criterion (currently disabled)
!   	if (particle(iCluster)%nPpxP > 1d0 .and. &
!   		.not. IsParticleStable(particle, nClusters, iCluster, Rho, waterDynVisco, kolmogorovLengthScale)) return

	! If we reached here, the cluster is OK for collision
  	eligible = .true.
  	
end function IsClusterIndividuallyEligibleForCollision

! ========================================================================================

logical function ShouldSkipCollision(particle, nClusters, iClusterA, iClusterB, &
	Rho, waterDynVisco, kolmogorovLengthScale) result(skip)

	!-------------------------------------------------------------------------------------	
	! Helper function for particle pair evaluation used in CoagulateParticles and CollideParticles.
	! See notes in above function IsClusterIndividuallyEligibleForCollision
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iClusterA, iClusterB
	real*8, intent(in) :: Rho, waterDynVisco, kolmogorovLengthScale
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle
	
    skip = .true.

	! 1) Require both clusters to be in phase 1. If either is not phase 1, skip.
    if (particle(iClusterA)%phase /= 1 .or. particle(iClusterB)%phase /= 1) return

	! 2) Require both clusters to have at least one particle 
    if (particle(iClusterA)%nPxC < 1d0 .or. particle(iClusterB)%nPxC < 1d0) return

	! 3) Require both clusters to meet operational size (here checking diameter >= min)
    if (2d0*particle(iClusterA)%radius < operational_size_poc_min) return
    if (2d0*particle(iClusterB)%radius < operational_size_poc_min) return

	! 4) Require both clusters to have nPpxP below the maximum allowed for collision
    if (MAX(particle(iClusterA)%nPpxP, particle(iClusterB)%nPpxP) >= maxPrimParticlesPerAggregateForCollision) return

!   Optional stability criterion (currently disabled)
! 	if (particle(iClusterA)%nPpxP > 1d0 .and. &
! 		.not. IsParticleStable(particle, nClusters, iClusterA, Rho, waterDynVisco, kolmogorovLengthScale)) return
!     if (particle(iClusterB)%nPpxP > 1d0 .and. &
!     	.not. IsParticleStable(particle, nClusters, iClusterB, Rho, waterDynVisco, kolmogorovLengthScale)) return
	
	! If we reached here, both clusters are OK for collision
  	skip = .false.

end function ShouldSkipCollision

! ========================================================================================

logical function IsParticleStable(particle, nClusters, iCluster, Rho, waterDynVisco, &
	kolmogorovLengthScale) result(stable)

	integer, intent(in) :: nClusters, iCluster
	real*8, intent(in) :: Rho, waterDynVisco, kolmogorovLengthScale
	type(lagrangianStateVars), dimension(nClusters), intent(in) :: particle

	real*8 :: ReFactor, ReParticle, ReTurb
	
    call CalculateParticleReynoldsNumber(ReFactor, ReParticle, particle, nClusters, iCluster, Rho, waterDynVisco)
	call CalculateTurbulenceScaleReynoldsNumber(ReTurb, particle, nClusters, iCluster, kolmogorovLengthScale)
    
    stable = (ReParticle <= breaking_reynolds_threshold .and. &
              ReTurb     <= breaking_reynolds_threshold)
                                  
end function IsParticleStable

! ========================================================================================

subroutine CollideParticles(particle, nClusters, iClusterA, iClusterB, nTotalPairs, &
	nPairsSubset, SMSterm, auxTerm, auxCount, TempK, Rho, waterDynVisco, waterShearRate, &
	kolmogorovLengthScale, gridCellVolume, iProfile)

	integer, intent(in) :: nClusters, iClusterA, iClusterB, nTotalPairs, nPairsSubset, iProfile
	real*8, intent(in) :: TempK, Rho, waterDynVisco, waterShearRate, kolmogorovLengthScale, &
		gridCellVolume
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm	
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	
	real*8 :: collisionTimeFrame, semiTime, radA, radB, radAB, veloA, veloB, mu, &
		collisionParam, brownianMotionKernel, fluidShearKernel, diffSettlingKernel, &
		collisionKernel, collisionEfficiency, coagulationProbability, harvest, &
		nParticlesInClusterWithMoreParticles, nParticlesInClusterWithLessParticles, &
		ejSampled, ratioMoreToLessParticles, ejRealised, ejExpected, remainder

	nTimesEnteringCollLoop(iProfile) = nTimesEnteringCollLoop(iProfile) + 1	
	auxCount(iCoagulationSuccess) = auxCount(iCoagulationSuccess) + particle(iClusterA)%nPxC + particle(iClusterB)%nPxC
	
	! ------ Initialise this ------
	collisionTimeFrame = timeStep ! s
	semiTime = 0 ! s
		
	do while (semiTime < timeStep) 		
		if (ShouldSkipCollision(particle, nClusters, iClusterA, iClusterB, Rho, waterDynVisco, kolmogorovLengthScale)) return ! exit this subroutine

		! semiTime is a timer and may become shorter than the loop time step when the
		! computed coagulation probability exceeds 1. In that case, the collision time
		! frame is reduced to scale the probability down. Coagulation is then evaluated
		! multiple times within this shorter time frame so that, when scaled back to the
		! original time step, the correct number of coagulation events is accounted for.
	
		! ------ Some unit adjustments ------
		radA = particle(iClusterA)%radius*1d-6 ! m
		radB = particle(iClusterB)%radius*1d-6 ! m
		radAB = radA + radB ! m
		veloA = particle(iClusterA)%velocity/SECONDS_PER_DAY ! m s-1
		veloB = particle(iClusterB)%velocity/SECONDS_PER_DAY ! m s-1
		mu = waterDynVisco/10d0 ! kg m-1 s-1
		
		if (.not. ieee_is_finite(radA) .or. .not. ieee_is_finite(radB) .or. radA <= 0d0 .or. radB <= 0d0) then
    		write(*,*) 'ERROR: invalid particle radii in CollideParticles'
    		write(*,*) '  radA, radB:', radA, radB
    		call WriteStatusAndStop()
		end if
			
		! ------- Build the collision kernel (m3 s-1) --------
		brownianMotionKernel = 2d0*BOLTZMANN_CNT*TempK*(radAB**2)/(3d0*mu*radA*radB) ! order 1e-17
		collisionParam = MIN(radA,radB) / MAX(radA,radB)        

		select case (choiceCollisionKernels)
		case (1) ! Burd & Jackson (2009)
	
			if (kolmogorovLengthScale > MAX(radA,radB)) then ! laminar flow
				fluidShearKernel = (4d0/3d0) * waterShearRate * radAB**3 		
			else ! turbulent flow
				select case (choiceIncludeHydrodynamicForces)
				case (1) ! include them
					fluidShearKernel = (9.8d0 * (collisionParam**2/(1d0 + 2d0*collisionParam**2)) * waterShearRate * radAB**3)
				case (0) ! do not
					fluidShearKernel = (1.3d0 * waterShearRate * radAB**3) 
				end select
			end if
			
			select case (choiceIncludeHydrodynamicForces)
			case (1) ! include them
				diffSettlingKernel = 0.5d0*PI*(MIN(radA,radB)**2)*ABS(veloA-veloB)	
			case (0)
				diffSettlingKernel = PI*(radAB**2)*ABS(veloA-veloB)
			end select
	
		case (2) ! Jackson (2001), with correction factors	
	
			if (kolmogorovLengthScale > MAX(radA,radB)) then ! laminar flow
				fluidShearKernel = (4d0/3d0) * waterShearRate * (radAB**3) * collisionParam**0.785d0
			else ! turbulent flow
				fluidShearKernel = 1.3d0 * waterShearRate * (radAB**3) * collisionParam**0.785d0
			end if
			
			if ((veloA-veloB)*(radA-radB) >= 0d0) then
				diffSettlingKernel = PI * (radAB**2) * ABS(veloA-veloB) * collisionParam**0.984d0
			else
				diffSettlingKernel = PI * (radAB**2) * ABS(veloA-veloB)
			end if
		
		end select		

		if (.not. ieee_is_finite(brownianMotionKernel) .or. &
    		.not. ieee_is_finite(fluidShearKernel) .or. &
    		.not. ieee_is_finite(diffSettlingKernel)) then
    		write(*,*) 'ERROR: non-finite collision kernel component'
    		write(*,*) '  Brownian, Shear, Settling:', brownianMotionKernel, fluidShearKernel, diffSettlingKernel
    		call WriteStatusAndStop()
		end if
	
		collisionKernel = brownianMotionKernel + fluidShearKernel + diffSettlingKernel ! m3 s-1
		
		! ------ Store diagnostics ------
		auxTerm(iCollisionKernel) = auxTerm(iCollisionKernel) + collisionKernel ! m3 s-1
		auxCount(iCollisionKernel) = auxCount(iCollisionKernel) + 1
		auxTerm(iBrownianKernel) = auxTerm(iBrownianKernel) + brownianMotionKernel ! m3 s-1
		auxCount(iBrownianKernel) = auxCount(iBrownianKernel) + 1
		auxTerm(iShearKernel) = auxTerm(iShearKernel) + fluidShearKernel ! m3 s-1
		auxCount(iShearKernel) = auxCount(iShearKernel) + 1
		auxTerm(iSettlingKernel) = auxTerm(iSettlingKernel) + diffSettlingKernel ! m3 s-1
		auxCount(iSettlingKernel) = auxCount(iSettlingKernel) + 1

		! -------- Collision efficiency, 0-1 (dimensionless) --------   
		collisionEfficiency = MAX(particle(iClusterA)%stickiness,particle(iClusterB)%stickiness)
		  
		! -------- Coagulation probability (dimensionless) --------
		! Scale coagulation probability by (nTotalPairs/nPairsSubset) to account for subsampling 
		! nPairsSubset out of nTotalPairs (in case subsampling was done)
		coagulationProbability = (dble(nTotalPairs) / dble(nPairsSubset)) &
			* collisionEfficiency * collisionKernel * collisionTimeFrame/gridCellVolume ! 0-1
		if (.not. ieee_is_finite(coagulationProbability)) then
    		write(*,*) 'ERROR: coagulationProbability not finite'
    		call WriteStatusAndStop()
		end if
		
		! NOTICE: the collision rate is very small (order 1e-14), the collision efficiency 
		! is between 0-1, the collision time frame is in the order 10-1000, and the gridCellVolume 
		! is 1 m3. Thus, the coagulation probability will be very small, <<1, but above 0 if 
		! at least one of the particles is sticky at all.	

		do while (coagulationProbability > 1d0)
			collisionTimeFrame = collisionTimeFrame * decreaseCollisionTimeframeFactor ! s
			coagulationProbability = (dble(nTotalPairs) / dble(nPairsSubset)) &
				* collisionEfficiency * collisionKernel * collisionTimeFrame/gridCellVolume
		end do
		semiTime = semiTime + collisionTimeFrame ! update timer
		
		! ------ Store diagnostics ------
		auxTerm(iCoagulationProbability) = auxTerm(iCoagulationProbability) + coagulationProbability
		auxCount(iCoagulationProbability) = auxCount(iCoagulationProbability) + 1
	
		! ------- Evaluate sticking --------
		if (coagulationProbability > 0d0) then
			nParticlesInClusterWithMoreParticles = MAX(particle(iClusterA)%nPxC, particle(iClusterB)%nPxC)
			nParticlesInClusterWithLessParticles = MIN(particle(iClusterA)%nPxC, particle(iClusterB)%nPxC)
			
			! Guard: counts must be > 0
			if (nParticlesInClusterWithMoreParticles <= 0d0 .or. nParticlesInClusterWithLessParticles <= 0d0) then
				write(*,*) 'ERROR: zero/invalid cluster count before ratio in CollideParticles.'
				write(*,*) '  semiTime:', semiTime
				write(*,*) '  Coagulation probability:', coagulationProbability
				write(*,*) '  nPxC (A,B):', particle(iClusterA)%nPxC, particle(iClusterB)%nPxC
				write(*,*) '  nPpxP (A,B):', particle(iClusterA)%nPpxP, particle(iClusterB)%nPpxP
				write(*,*) '  Depth (A,B):', particle(iClusterA)%depth, particle(iClusterB)%depth
				write(*,*) '  Porosity (A,B):', particle(iClusterA)%porosity, particle(iClusterB)%porosity
				write(*,*) '  Mass (A,B):', particle(iClusterA)%mass, particle(iClusterB)%mass
				write(*,*) '  Radius (A,B):', particle(iClusterA)%radius, particle(iClusterB)%radius
				write(*,*) '  Velocity (A,B):', particle(iClusterA)%velocity, particle(iClusterB)%velocity
				write(*,*) '  Particle type (A,B):', particle(iClusterA)%initType, particle(iClusterB)%initType
				write(*,*) '  Particle phase (A,B):', particle(iClusterA)%phase, particle(iClusterB)%phase
				call WriteStatusAndStop()
			end if

			auxTerm(iNumPxCmore) = auxTerm(iNumPxCmore) + nParticlesInClusterWithMoreParticles
			auxCount(iNumPxCmore) = auxCount(iNumPxCmore) + 1
			auxTerm(iNumPxCless) = auxTerm(iNumPxCless) + nParticlesInClusterWithLessParticles
			auxCount(iNumPxCless) = auxCount(iNumPxCless) + 1
			
			! ------ Calculate the number of particles exchanged ------
			! Ej is the no. particles of one cluster that collide with one particle of 
			! the other cluster. Currently, it is a combinatorial expectation derived from considering 
			! all pairwise possibilities and how many of the n_i particles collide with any 
			! of the n_j particles
			ratioMoreToLessParticles = nParticlesInClusterWithMoreParticles / nParticlesInClusterWithLessParticles
			if (coagulationProbability >= 1d0) then ! --> everything collides
    			ejExpected = ratioMoreToLessParticles
			else if (coagulationProbability <= 0d0) then ! --> nothing collides
    			ejExpected = 0d0
			else
    			ejExpected = ratioMoreToLessParticles * &
             		(1d0 - (1d0 - coagulationProbability)**nParticlesInClusterWithLessParticles)
			end if
			if (.not. ieee_is_finite(ejExpected) .or. ejExpected < 0d0) then
    			write(*,*) 'ERROR: invalid ejSampled in CollideParticles'
    			write(*,*) '  ejExpected:', ejExpected
    			write(*,*) '  ratioMoreToLessParticles:', ratioMoreToLessParticles
    			write(*,*) '  coagulationProbability:', coagulationProbability
    			call WriteStatusAndStop()
			end if

			! ------ Stochastic rounding ------
			! We need to make the number of particles an integer number. We round the number of 
			! particles up or down by comparing the fractional part of the number with a random number
			ejRealised = SafeFloorNonNegative(ejExpected) ! round down
			remainder = ejExpected - ejRealised ! used as the probability of adding one extra particle
			remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
			call random_number(harvest)
			if (harvest < remainder) then
				ejSampled = ejRealised + 1d0
			else
				ejSampled = ejRealised
			end if
			if (.not. ieee_is_finite(ejSampled) .or. ejSampled < 0d0) then
    			write(*,*) 'ERROR: invalid ejSampled after stochastic rounding'
    			write(*,*) '  ejExpected:', ejExpected
    			write(*,*) '  ejSampled:', ejSampled
    			call WriteStatusAndStop()
			end if
					  
			! ------ Realise particle sticking ------
			if (ejSampled >= 1d0) then
				call StickParticles(particle, nClusters, iClusterA, iClusterB, SMSterm, auxTerm, &
					ejSampled, Rho, waterDynVisco, iProfile)
			end if
	
		end if ! coagulationProbability > 0
	end do ! semiTime while loop 

end subroutine CollideParticles

! ========================================================================================

subroutine StickParticles(particle, nClusters, iClusterA, iClusterB, SMSterm, auxTerm, &
	ejSampled, Rho, waterDynVisco, iProfile)  

	!-------------------------------------------------------------------------------------	
	! There are 4 collision types.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, iClusterA, iClusterB, iProfile
	real*8, intent(in) :: ejSampled, Rho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm
	
	integer :: colliType
	real*8 :: nParticlesPerClusterA, nParticlesPerClusterB, molesOrgCbeforeAggregation, &
		molesOrgCafterAggregation, radAbefore, radBbefore, depthAbefore, depthBbefore, &
		nPrimParticlesPerParticleA, nPrimParticlesPerParticleB, totalParticles, nA, nB, ejClamped
	real*8, parameter :: EQ_EPS = 1d-12
				
	nParticlesPerClusterA = particle(iClusterA)%nPxC
	nParticlesPerClusterB = particle(iClusterB)%nPxC	
	molesOrgCbeforeAggregation = particle(iClusterA)%molesOrgC*particle(iClusterA)%nPxC &
		+ particle(iClusterB)%molesOrgC*particle(iClusterB)%nPxC
	radAbefore = particle(iClusterA)%radius
	radBbefore = particle(iClusterB)%radius
	depthAbefore = particle(iClusterA)%depth
	depthBbefore = particle(iClusterB)%depth
	nPrimParticlesPerParticleA = particle(iClusterA)%nPpxP
	nPrimParticlesPerParticleB = particle(iClusterB)%nPpxP
				
	auxTerm(iCoagulationSuccess) = auxTerm(iCoagulationSuccess) + nParticlesPerClusterA + nParticlesPerClusterB
	
	colliType = 0
	
	!-------------------------------------------------------------------------------------
	! Type 1 collision
	! Cluster A looses particles, which belong now to cluster B. This implies that 
	! cluster A looses particles per cluster, but remains with the same number of primary 
	! particles per particle, and cluster B gains primary particles per particle, but 
	! remains with the same number of particles per cluster. Only the biogeochemical 
	! attributes of cluster B change.
	!-------------------------------------------------------------------------------------
	
	if (nParticlesPerClusterA > nParticlesPerClusterB .and. &
		(ejSampled*nParticlesPerClusterB) < nParticlesPerClusterA) then
	
		! Ensure donor keeps at least one particle; if not, abort sticking completely
		ejClamped = MaxIntegerExchangeFactor(ejSampled, nParticlesPerClusterA, nParticlesPerClusterB)
		if (ejClamped <= 0d0) return
	
		nCollisionsType1(iProfile) = nCollisionsType1(iProfile) + 1

		particle(iClusterA)%nPxC = particle(iClusterA)%nPxC - (ejClamped * particle(iClusterB)%nPxC)
		particle(iClusterB)%nPpxP = particle(iClusterB)%nPpxP + (ejClamped * particle(iClusterA)%nPpxP)
		
		if (particle(iClusterA)%nPxC < 1d0) then
    		write(*,*) 'ERROR: donor cluster emptied in collision type 1'
    		write(*,*) '  nPxC_B =', particle(iClusterA)%nPxC
    		write(*,*) '  ej =', ejClamped
    		call WriteStatusAndStop()
		end if

		particle(iClusterA)%nPpxC = particle(iClusterA)%nPxC * particle(iClusterA)%nPpxP
		particle(iClusterB)%nPpxC = particle(iClusterB)%nPxC * particle(iClusterB)%nPpxP

		particle(iClusterB)%molesOrgC = particle(iClusterB)%molesOrgC + (ejClamped * particle(iClusterA)%molesOrgC)
		particle(iClusterB)%molesTepC = particle(iClusterB)%molesTepC + (ejClamped * particle(iClusterA)%molesTepC)
		particle(iClusterB)%molesMineral(:) = particle(iClusterB)%molesMineral(:) + (ejClamped * particle(iClusterA)%molesMineral(:))
	
		! Sanity check
		molesOrgCafterAggregation = particle(iClusterA)%molesOrgC*particle(iClusterA)%nPxC &
								  + (particle(iClusterB)%molesOrgC*particle(iClusterB)%nPxC)
								  
		call PostParticleDynamicsMassConservationCheck(SMSterm, molesOrgCafterAggregation, molesOrgCbeforeAggregation, &
			nParticlesPerClusterA, nParticlesPerClusterB, radAbefore, radBbefore, 'collision type 1')
		
		! Update cluster attributes	
		particle(iClusterB)%faecal = 0
		particle(iClusterB)%living = 0
		particle(iClusterB)%initType = ParticleUpdatedType(particle, nClusters, iClusterB)
		call ParticleFractalDimension(particle, nClusters, iClusterB, particle(iClusterB)%initType)
		call ParticleDryMass(particle, nClusters, iClusterB)
		call ParticleMaterialVolume(particle, nClusters, iClusterB)
		call ParticleRadius(particle, nClusters, iClusterB)
		call ParticlePorosity(particle, nClusters, iClusterB)
		call ParticleDensity(particle, nClusters, iClusterB, Rho)
		call ParticleStickiness(particle, nClusters, iClusterB)	
		call ParticleSettlingVelocity(particle, nClusters, iClusterB, Rho, waterDynVisco)
		
		! Check: is the particle physically consistent?
		call CheckParticleSanity(particle, nClusters, iClusterA, 'status A after collision') 
		call CheckParticleSanity(particle, nClusters, iClusterB, 'status B after collision')

	!-------------------------------------------------------------------------------------
	! Type 2 collision
	! Cluster B looses particles, which belong now to cluster A. This implies that 
	! cluster B looses particles per cluster, but remains with the same number of primary 
	! particles per particle, and cluster A gains primary particles per particle, but 
	! remains with the same number of particles per cluster. Only the biogeochemical 
	! attributes of cluster A change.
	!-------------------------------------------------------------------------------------
	
	else if (nParticlesPerClusterA < nParticlesPerClusterB .and. &
		(ejSampled*nParticlesPerClusterA) < nParticlesPerClusterB) then

		! Ensure donor keeps at least one particle; if not, abort sticking completely
		ejClamped = MaxIntegerExchangeFactor(ejSampled, nParticlesPerClusterB, nParticlesPerClusterA)		
		if (ejClamped <= 0d0) return
		
		colliType = 2
		nCollisionsType2(iProfile) = nCollisionsType2(iProfile) + 1

		particle(iClusterB)%nPxC = particle(iClusterB)%nPxC - (ejClamped * particle(iClusterA)%nPxC)
		particle(iClusterA)%nPpxP = particle(iClusterA)%nPpxP + (ejClamped * particle(iClusterB)%nPpxP)
		
		if (particle(iClusterB)%nPxC < 1d0) then
    		write(*,*) 'ERROR: donor cluster emptied in collision type 2'
    		write(*,*) '  nPxC_B =', particle(iClusterB)%nPxC
    		write(*,*) '  ejClamped =', ejClamped
    		call WriteStatusAndStop()
		end if

		particle(iClusterA)%nPpxC = particle(iClusterA)%nPxC * particle(iClusterA)%nPpxP
		particle(iClusterB)%nPpxC = particle(iClusterB)%nPxC * particle(iClusterB)%nPpxP
			
		particle(iClusterA)%molesOrgC = particle(iClusterA)%molesOrgC + (ejClamped * particle(iClusterB)%molesOrgC)
		particle(iClusterA)%molesTepC = particle(iClusterA)%molesTepC + (ejClamped * particle(iClusterB)%molesTepC)
		particle(iClusterA)%molesMineral(:) = particle(iClusterA)%molesMineral(:) + (ejClamped * particle(iClusterB)%molesMineral(:))

		! Sanity check
		molesOrgCafterAggregation = particle(iClusterA)%molesOrgC*particle(iClusterA)%nPxC &
								  + particle(iClusterB)%molesOrgC*particle(iClusterB)%nPxC
		call PostParticleDynamicsMassConservationCheck(SMSterm, molesOrgCafterAggregation, molesOrgCbeforeAggregation, &
			nParticlesPerClusterA, nParticlesPerClusterB, radAbefore, radBbefore, 'collision type 2')
		
		! Update cluster A attributes	
		particle(iClusterA)%faecal = 0		
		particle(iClusterA)%living = 0
		particle(iClusterA)%initType = ParticleUpdatedType(particle, nClusters, iClusterA)
		call ParticleFractalDimension(particle, nClusters, iClusterA, particle(iClusterA)%initType)
		call ParticleDryMass(particle, nClusters, iClusterA) 
		call ParticleMaterialVolume(particle, nClusters, iClusterA) 
		call ParticleRadius(particle, nClusters, iClusterA)
		call ParticlePorosity(particle, nClusters, iClusterA)
		call ParticleDensity(particle, nClusters, iClusterA, Rho) 
		call ParticleStickiness(particle, nClusters, iClusterA) 
		call ParticleSettlingVelocity(particle, nClusters, iClusterA, Rho, waterDynVisco)
			
		! Check: is the particle physically consistent?
		call CheckParticleSanity(particle, nClusters, iClusterA, 'status A after collision') 
		call CheckParticleSanity(particle, nClusters, iClusterB, 'status B after collision')
		
	!-------------------------------------------------------------------------------------
	! Type 3 collision
	! We create a new cluster A and a new cluster B. Each particle of cluster A sticks to
	! each particle of cluster B. We then split the new aggregates created into two new 
	! clusters, the new cluster A and the new cluster B. The biogeochemical attributes of 
	! both clusters are updated.
	!-------------------------------------------------------------------------------------

	else if (ABS(nParticlesPerClusterA - nParticlesPerClusterB) <= EQ_EPS .and. nParticlesPerClusterA > 1d0) then

		colliType = 3
		nCollisionsType3(iProfile) = nCollisionsType3(iProfile) + 1
	
		particle(iClusterA)%nPpxP = particle(iClusterA)%nPpxP + particle(iClusterB)%nPpxP
		particle(iClusterB)%nPpxP = particle(iClusterA)%nPpxP ! =nPxP count of A
		
		! Split particles
		totalParticles = nParticlesPerClusterA
		call SafeSplitIntegerLike(totalParticles, nA, nB)
		particle(iClusterA)%nPxC = nA
		particle(iClusterB)%nPxC = nB

		particle(iClusterA)%nPpxC = particle(iClusterA)%nPxC * particle(iClusterA)%nPpxP
		particle(iClusterB)%nPpxC = particle(iClusterB)%nPxC * particle(iClusterB)%nPpxP

		particle(iClusterA)%molesOrgC = particle(iClusterA)%molesOrgC + particle(iClusterB)%molesOrgC
		particle(iClusterB)%molesOrgC = particle(iClusterA)%molesOrgC
		particle(iClusterA)%molesTepC = particle(iClusterA)%molesTepC + particle(iClusterB)%molesTepC
		particle(iClusterB)%molesTepC = particle(iClusterA)%molesTepC
		particle(iClusterA)%molesMineral(:) = particle(iClusterA)%molesMineral(:) + particle(iClusterB)%molesMineral(:)
		particle(iClusterB)%molesMineral(:) = particle(iClusterA)%molesMineral(:)

		! Sanity check
		molesOrgCafterAggregation = particle(iClusterA)%molesOrgC*particle(iClusterA)%nPxC &
								  + particle(iClusterB)%molesOrgC*particle(iClusterB)%nPxC
		call PostParticleDynamicsMassConservationCheck(SMSterm, molesOrgCafterAggregation, molesOrgCbeforeAggregation, &
			nParticlesPerClusterA, nParticlesPerClusterB, radAbefore, radBbefore, 'collision type 3')
		
		! Update both cluster attributes
		particle(iClusterA)%faecal = 0		
		particle(iClusterA)%living = 0
		particle(iClusterA)%initType = ParticleUpdatedType(particle, nClusters, iClusterA)
		call ParticleFractalDimension(particle, nClusters, iClusterA, particle(iClusterA)%initType)
		call ParticleDryMass(particle, nClusters, iClusterA)
		call ParticleMaterialVolume(particle, nClusters, iClusterA)
		call ParticleRadius(particle, nClusters, iClusterA)
		call ParticlePorosity(particle, nClusters, iClusterA)
		call ParticleDensity(particle, nClusters, iClusterA, Rho)
		call ParticleStickiness(particle, nClusters, iClusterA) 	
		call ParticleSettlingVelocity(particle, nClusters, iClusterA, Rho, waterDynVisco)
		
		particle(iClusterB)%faecal = 0
		particle(iClusterB)%living = 0 
		particle(iClusterB)%initType = ParticleUpdatedType(particle, nClusters, iClusterB)
		call ParticleFractalDimension(particle, nClusters, iClusterB, particle(iClusterB)%initType)
		call ParticleDryMass(particle, nClusters, iClusterB)
		call ParticleMaterialVolume(particle, nClusters, iClusterB)
		call ParticleRadius(particle, nClusters, iClusterB)
		call ParticlePorosity(particle, nClusters, iClusterB)
		call ParticleDensity(particle, nClusters, iClusterB, Rho)
		call ParticleStickiness(particle, nClusters, iClusterB) 
		call ParticleSettlingVelocity(particle, nClusters, iClusterB, Rho, waterDynVisco)
		
		! Check: is the particle physically consistent?
		call CheckParticleSanity(particle, nClusters, iClusterA, 'status A after collision') 
		call CheckParticleSanity(particle, nClusters, iClusterB, 'status B after collision')
		
	!-------------------------------------------------------------------------------------
	! Type 4 collision
	! When both clusters are formed just by one particle only, we coagulate both particles
	! in one cluster (A) and empty the other cluster (B), which takes phase 3.
	!-------------------------------------------------------------------------------------

	else if (nParticlesPerClusterA == nParticlesPerClusterB .and. nParticlesPerClusterA == 1d0) then

		colliType = 4
		nCollisionsType4(iProfile) = nCollisionsType4(iProfile) + 1
		
		particle(iClusterA)%nPpxP = particle(iClusterA)%nPpxP + particle(iClusterB)%nPpxP
		particle(iClusterA)%nPpxC = particle(iClusterA)%nPxC * particle(iClusterA)%nPpxP

		particle(iClusterA)%molesOrgC = particle(iClusterA)%molesOrgC + particle(iClusterB)%molesOrgC
		particle(iClusterA)%molesTepC = particle(iClusterA)%molesTepC + particle(iClusterB)%molesTepC
		particle(iClusterA)%molesMineral(:) = particle(iClusterA)%molesMineral(:) + particle(iClusterB)%molesMineral(:)

		! Sanity check
		molesOrgCafterAggregation = particle(iClusterA)%molesOrgC*particle(iClusterA)%nPxC
		call PostParticleDynamicsMassConservationCheck(SMSterm, molesOrgCafterAggregation, molesOrgCbeforeAggregation, &
			nParticlesPerClusterA, nParticlesPerClusterB, radAbefore, radBbefore, 'collision type 4')
		
		! Update cluster A attributes
		particle(iClusterA)%faecal = 0		
		particle(iClusterA)%living = 0
		particle(iClusterA)%initType = ParticleUpdatedType(particle, nClusters, iClusterA)
		call ParticleFractalDimension(particle, nClusters, iClusterA, particle(iClusterA)%initType)
		call ParticleDryMass(particle, nClusters, iClusterA)
		call ParticleMaterialVolume(particle, nClusters, iClusterA)
		call ParticleRadius(particle, nClusters, iClusterA)
		call ParticlePorosity(particle, nClusters, iClusterA)
		call ParticleDensity(particle, nClusters, iClusterA, Rho)
		call ParticleStickiness(particle, nClusters, iClusterA) 
		call ParticleSettlingVelocity(particle, nClusters, iClusterA, Rho, waterDynVisco)			

		! Send cluster B to empty phase
		particle(iClusterB)%phase = 3
		particle(iClusterB)%nPxC = 0d0 
		
		! Check: is the particle physically consistent?
		call CheckParticleSanity(particle, nClusters, iClusterA, 'status A after collision') 

	end if

end subroutine StickParticles

! ========================================================================================

function MaxIntegerExchangeFactor(nDonorParticlesThatCollideWithOneReceiverParticle, &
	nDonorParticlesPerCluster, nReceiverParticlesPerCluster) result(ej)
	
	!-------------------------------------------------------------------------------------
	! This function ensures a donor cluster must never be emptied, i.e., it must retain at 
	! least one particle per cluster. It returns maximum integer number of 
	! donor-particles-per-receiver-particle that can be exchanged without emptying the 
	! donor cluster.
  	!-------------------------------------------------------------------------------------
  	
  	real*8, intent(in) :: nDonorParticlesThatCollideWithOneReceiverParticle, &
		nDonorParticlesPerCluster, nReceiverParticlesPerCluster
  	real*8 :: ej, ejMax

	!  Sanitise inputs
  	if (nReceiverParticlesPerCluster <= 0d0) then
    	ej = 0d0
    	return
  	end if
  	if (nDonorParticlesPerCluster <= 1d0) then
  		ej = 0d0
  		return 
  	end if
  	if (nDonorParticlesThatCollideWithOneReceiverParticle <= 0d0) then
  		ej = 0d0
  		return 
  	end if
  	
  	! NaN and Inf checks
    if ( nDonorParticlesThatCollideWithOneReceiverParticle /= &
         nDonorParticlesThatCollideWithOneReceiverParticle .or. &
         nDonorParticlesPerCluster /= nDonorParticlesPerCluster .or. &
         nReceiverParticlesPerCluster /= nReceiverParticlesPerCluster ) then
        write(*,*) 'ERROR: NaN in MaxIntegerExchangeFactor'
        call WriteStatusAndStop()
    end if
    if ( abs(nDonorParticlesThatCollideWithOneReceiverParticle) > huge(1d0) .or. &
     	abs(nDonorParticlesPerCluster) > huge(1d0) .or. &
     	abs(nReceiverParticlesPerCluster) > huge(1d0) ) then
    	write(*,*) 'ERROR: Inf in MaxIntegerExchangeFactor'
    	call WriteStatusAndStop()
	end if
    
    ! -------- Compute Ej_max safely --------
    ejMax = (nDonorParticlesPerCluster - 1d0) / nReceiverParticlesPerCluster

    if (ejMax <= 0d0) then
        ej = 0d0
        return
    end if

    ! REAL floor, no INTEGER conversion
    ejMax = AINT(ejMax)

    ! -------- Apply bound --------
    ej = MIN(nDonorParticlesThatCollideWithOneReceiverParticle, ejMax)
    
    if (abs(ej - AINT(ej)) > 1d-10) then
    	write(*,*) 'ERROR: ej not integer-like'
    	write(*,*) '  ej =', ej
    	call WriteStatusAndStop()
	end if
  
end function MaxIntegerExchangeFactor

! ========================================================================================

subroutine BreakUnstableParticles(particle, nClusters, nLayerClusters, layerClusterIndices, &
	SMSterm, Rho, waterDynVisco, kolmogorovLengthScale, iProfile)

	integer, intent(in) :: nClusters, nLayerClusters, iProfile
	integer, dimension(nLayerClusters), intent(in) :: layerClusterIndices	
	real*8, intent(in) :: Rho, waterDynVisco, kolmogorovLengthScale
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm	
	
	integer :: iLayerCluster, iCluster, nBreaks, maxBreaks
	real*8 :: nPrimParticles, radiusP, log2n
	logical :: shouldBreak

	do iLayerCluster = 1, nLayerClusters	
		iCluster = layerClusterIndices(iLayerCluster)
		nPrimParticles = particle(iCluster)%nPpxP
        radiusP = particle(iCluster)%radius
		
    	if (particle(iCluster)%phase /= 1) cycle
    	if (nPrimParticles < 2d0) cycle
    	select case (choiceCriteriaToBreak)
    	case (1)
			if (radiusP <= kolmogorovLengthScale*1d6) cycle
		case (2)
			if (IsParticleStable(particle, nClusters, iCluster, Rho, waterDynVisco, kolmogorovLengthScale)) cycle
		end select
		
		if (.not. ieee_is_finite(nPrimParticles)) then
    		write(*,*) 'ERROR: invalid nPpxP before fragmentation:', nPrimParticles
    		call WriteStatusAndStop()
		end if
		
		! -------- Determine a safe upper bound on number of fragmentation events --------
    	! Each fragmentation event is assumed to reduce aggregate size by at least a factor 
    	! of ~2 (i.e., each break produces ≥2 daughter particles). At most log2(nPpxP) 
    	! fragmentation events are needed to reach monomers. Add +3 as a conservative safety margin.
    	! “How many times can I divide nPpxP by 2 until I reach ~1?
    	log2n = SafeLog(nPrimParticles) / log(2d0)
    	if (.not. ieee_is_finite(log2n) .or. log2n < 0d0) then
    		write(*,*) 'ERROR: invalid log2(nPpxP):', log2n
    		call WriteStatusAndStop()
		end if
  		
  		maxBreaks = SafeCeilingNonNegative(log2n) + 3
		nBreaks = 0
		
		! ---------- Fragment until stable or limit reached ----------
		do
			! Stop if cluster no longer eligible
			if (particle(iCluster)%phase /= 1) exit
			if (particle(iCluster)%nPpxP < 2d0) exit
			select case (choiceCriteriaToBreak)
			case (1)
				shouldBreak = (particle(iCluster)%radius > kolmogorovLengthScale * 1d6)
			case (2)
				shouldBreak = .not. IsParticleStable(particle, nClusters, iCluster, Rho, waterDynVisco, kolmogorovLengthScale)
			end select
		
			if (.not. shouldBreak) exit
			if (nBreaks >= maxBreaks) exit ! prevent infinite loop
				
			call BreakParticle(particle, nClusters, iCluster, SMSterm, Rho, waterDynVisco)
			
			nFragmentedBigClusters(iProfile) = nFragmentedBigClusters(iProfile) + 1
			nBreaks = nBreaks + 1
		end do

	end do

end subroutine BreakUnstableParticles

! ========================================================================================

subroutine BreakParticle(particle, nClusters, iCluster, SMSterm, Rho, waterDynVisco)

	!-------------------------------------------------------------------------------------	
	! The number of daughter particles is the number into which a particle of the cluster 
	! fragments. The following keeps nDaughters between 2–100 (as suggested by Nathan Briggs).
	! 100 is scaled down by the particle fragility, which is a function of the radius of the 
	! particle (larger particles are more fragile) and the amount of organic carbon (it keeps 
	! a particle more cohesive). We need to get real*8 variables with enforced discreteness 
	! via stochastic rounding
	!-------------------------------------------------------------------------------------	
  	
  	integer, intent(in) :: nClusters, iCluster
  	real*8, intent(in) :: Rho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm		

	real*8 :: nDaughters, nPpxParentParticle, nPxParentCluster, nPpxDaughterParticle, &
		nPxDaughterCluster, radiusParentParticle, fracOrgMatter, particleFragility, &
		parentClusterMolesOrgC, parentClusterMolesTepC, daughterClusterMolesOrgC, &
		maxNumDaughters, expectedNumDaughters, realisedNumDaughters, remainder, harvest, &
		expectedNumPrimPartPerDaughter, realisedNumPrimPartPerDaughter
	real*8, dimension(nMinerals) :: parentClusterMolesMineral

	nPpxParentParticle = particle(iCluster)%nPpxP
  	nPxParentCluster  = particle(iCluster)%nPxC
  	if (nPpxParentParticle < 2d0) return ! early return for monomers
  	
  	! ------ Sanity checks first ------
  	if (.not. ieee_is_finite(nPpxParentParticle) .or. nPpxParentParticle < 2d0) then
    	write(*,*) 'ERROR: invalid nPpxP before breaking', nPpxParentParticle
    	call WriteStatusAndStop()
	end if
  	if (.not. ieee_is_finite(nPxParentCluster) .or. nPxParentCluster <= 0d0) then
    	write(*,*) 'ERROR: invalid nPxC before breaking', nPxParentCluster
    	call WriteStatusAndStop()
	end if
	if (.not. ieee_is_finite(particle(iCluster)%radius) .or. particle(iCluster)%radius <= 0d0) then
    	write(*,*) 'ERROR: invalid radius before breaking', particle(iCluster)%radius
    	call WriteStatusAndStop()
	end if

	! ------ Compute the fragility factor ------
	radiusParentParticle = particle(iCluster)%radius
	fracOrgMatter = VolumetricFractionOfOrganicMatter( particle, nClusters, iCluster )
	particleFragility = (1d0 - fracOrgMatter) * MIN(radiusParentParticle,1d4)/1d4 ! 1=very fragile, 0=no fragile	
	maxNumDaughters = 100d0 * particleFragility
	
	! ------ Compute the number of daughters with stochastic rounding ------
	expectedNumDaughters = MIN(MAX(maxNumDaughters, 2d0), nPpxParentParticle) ! expected number of daughters (real)
	realisedNumDaughters = SafeFloorNonNegative(expectedNumDaughters) ! round down
	remainder = expectedNumDaughters - realisedNumDaughters ! used as the probability of adding one extra particle
	remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
	
	call random_number(harvest)
	if (harvest < remainder) then
		nDaughters = realisedNumDaughters + 1d0
	else
		nDaughters = realisedNumDaughters
	end if
	if (nDaughters < 2d0) then
		write(*,*) 'ERROR: nDaughters < 2 in BreakParticle'
		call WriteStatusAndStop()
	end if

	if (nDaughters >= 2d0 .and. nPpxParentParticle >= nDaughters) then

   		! ------ Calculate total primary particle count with stochastic rounding ------
		expectedNumPrimPartPerDaughter = SafeDivide(nPpxParentParticle, nDaughters, 0d0)
		realisedNumPrimPartPerDaughter = SafeFloorNonNegative(expectedNumPrimPartPerDaughter) ! round down
		remainder = expectedNumPrimPartPerDaughter - realisedNumPrimPartPerDaughter ! used as the probability of adding one extra particle
		remainder = MAX(0d0, MIN(1d0, remainder)) ! numerical safety
		
		call random_number(harvest)
		if (harvest < remainder) then
			nPpxDaughterParticle = realisedNumPrimPartPerDaughter + 1d0
		else
			nPpxDaughterParticle = realisedNumPrimPartPerDaughter
		end if
		if (nPpxDaughterParticle < 1d0) then
			write(*,*) 'ERROR: nPpxDaughterParticle < 1 in BreakParticle'
			call WriteStatusAndStop()
		end if

		! ------ Redistribution of particles in the daughter particle ------
		nPxDaughterCluster = nPxParentCluster * nDaughters
		if (.not. ieee_is_finite(nPxDaughterCluster) .or. nPxDaughterCluster <= 0d0) then
    		write(*,*) 'ERROR: invalid nPxDaughterCluster in BreakParticle', nPxDaughterCluster
    		call WriteStatusAndStop()
		end if
		
		particle(iCluster)%nPxC = nPxDaughterCluster
		particle(iCluster)%nPpxP = nPpxDaughterParticle
		particle(iCluster)%nPpxC = particle(iCluster)%nPxC * particle(iCluster)%nPpxP
	
		! ------ Redistribution of material in the daughter particle ------	
		parentClusterMolesOrgC = particle(iCluster)%molesOrgC * nPxParentCluster 				
		parentClusterMolesTepC = particle(iCluster)%molesTepC * nPxParentCluster 				
		parentClusterMolesMineral(:) = particle(iCluster)%molesMineral(:) * nPxParentCluster 
				
		particle(iCluster)%molesOrgC = parentClusterMolesOrgC / nPxDaughterCluster
		particle(iCluster)%molesTepC = parentClusterMolesTepC / nPxDaughterCluster
		particle(iCluster)%molesMineral(:) = parentClusterMolesMineral(:) / nPxDaughterCluster

		! ------ Update daughter attributes ------
		particle(iCluster)%living = 0
		!particle(iCluster)%faecal = 0 ! comment in --> if a faecal pellet breaks into two, it will keep being a fp
		particle(iCluster)%initType = ParticleUpdatedType(particle, nClusters, iCluster)	
		call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
		call ParticleDryMass(particle, nClusters, iCluster) 
		call ParticleMaterialVolume(particle, nClusters, iCluster) 
		call ParticleRadius(particle, nClusters, iCluster)
		call ParticlePorosity(particle, nClusters, iCluster)
		call ParticleDensity(particle, nClusters, iCluster, Rho)
		call ParticleStickiness(particle, nClusters, iCluster)	
		call ParticleSettlingVelocity(particle, nClusters, iCluster, Rho, waterDynVisco)

		! ------ Sanity checks ------
		daughterClusterMolesOrgC = particle(iCluster)%molesOrgC*particle(iCluster)%nPxC
		
		! Does the dauhter particle contain more material than the parent particle?
		if (IsExceedingInitialAmount(daughterClusterMolesOrgC, parentClusterMolesOrgC, detection_limit_poc)) then
			write(*,*) 'ERROR: more material is formed after breaking than initially available.'
			write(*,*) '  Starting and final material:', parentClusterMolesOrgC, daughterClusterMolesOrgC
			write(*,*) '  Radius:', particle(iCluster)%radius
			write(*,*) '  Mineral content:', particle(iCluster)%molesMineral(:)
			write(*,*) '  Depth:', particle(iCluster)%depth
			call WriteStatusAndStop()
		end if
		
		! If no mass conservation after breaking --> send leftovers to solubilisation pool
		if (.not. NearlyEqual(parentClusterMolesOrgC, daughterClusterMolesOrgC, detection_limit_poc)) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + (parentClusterMolesOrgC - daughterClusterMolesOrgC)
		end if
				
	end if
		
end subroutine BreakParticle

! ========================================================================================

subroutine PostParticleDynamicsMassConservationCheck(SMSterm, clusterOrgCmolAfter, &
	clusterOrgCmolBefore, nParticlesPerClusterA, nParticlesPerClusterB, radAbefore, &
	radBbefore, processDescription)
	
	! ------------------------------------------------------------------------------------
    ! Subroutine that supports StickParticles and PackSingleParticles.
    ! ------------------------------------------------------------------------------------
    
    real*8, intent(in) :: nParticlesPerClusterA, nParticlesPerClusterB, clusterOrgCmolAfter, &
    	clusterOrgCmolBefore, radAbefore, radBbefore
    character(len=*), intent(in), optional :: processDescription
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	
	character(len=:), allocatable :: procLabel
	
	if (present(processDescription)) then
      	procLabel = trim(processDescription)
    else
      	procLabel = 'unspecified'
    end if
	
	! Checking after =/ before
	if (IsExceedingInitialAmount(clusterOrgCmolAfter, clusterOrgCmolBefore, detection_limit_poc)) then
		write(*,*) 'ERROR: more material after ' // trim(procLabel) // ' than the starting one.'
		write(*,*) '  Before and after material:', clusterOrgCmolBefore, clusterOrgCmolAfter
		write(*,*) '  Before nPxC (A,B):', nParticlesPerClusterA, nParticlesPerClusterB
		write(*,*) '  Before radius (A,B):', radAbefore, radBbefore
		call WriteStatusAndStop()
	end if

	! Mass conservation sanity check --> redirect leftovers to solubilisation pool
	if (.not. NearlyEqual(clusterOrgCmolBefore, clusterOrgCmolAfter, detection_limit_poc)) then
		SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + (clusterOrgCmolBefore - clusterOrgCmolAfter)
	end if

end subroutine PostParticleDynamicsMassConservationCheck

! ========================================================================================

end module particledynamics