module mineraldissolution

! ----------------------------------------------------------------------------------------
! This module handles the dissolution of mineral components of particles (CaCO3 and opal)
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use modelcounters, only: nMineralClustersDissolved
use modelconstants, only: SECONDS_PER_DAY, MOLAR_MASS_OPAL, MOLAR_MASS_CACO3, RHO_OPAL, RHO_CALCITE
use modelparameters, only: timeStep, iCalcite, iOpal, iClay, iDissolCaCO3, iDissolOpal, &
	iMicrobSolubOrgC, iMicrobSolubTepC, iMicrobSolubClay, dissol_rate_calc, order_react_dissol_calc, &
	dissol_rate_opal_0deg, q10_bSi, maxNumSmsTerms, detection_limit_opal, detection_limit_calc, &
	detection_limit_clay, shrinkAfterMineralDissolutionScheme, operational_size_poc_min, detection_limit_poc
use sanitychecks, only: IsQuantityEffectivelyZero, IsExceedingInitialAmount, &
	WriteStatusAndStop, CheckParticleSanity
use calcparticleattributes, only: ParticleFractalDimension, ParticleDryMass, &
	ParticleMaterialVolume, ParticleRadius, ParticlePorosity, ParticleDensity, &
	ParticleStickiness, ParticleSettlingVelocity, ShrinkParticle
	
implicit none
private
public :: CalciumCarbonateAndOpalDissolution, CaCO3dissolutionLimFactor

contains

! ========================================================================================

subroutine CalciumCarbonateAndOpalDissolution(particle, nClusters, SMSterm, nLayerClusters, &
	layerClusterIndices, OmegaCalc, TempC, waterRho, waterDynVisco, iProfile)

	integer, intent(in) :: nClusters, nLayerClusters, iProfile	
	integer, dimension(nLayerClusters), intent(in) :: layerClusterIndices
	real*8, intent(in) :: OmegaCalc, TempC, waterRho, waterDynVisco
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	real*8, dimension(maxNumSmsTerms), intent(inout) :: SMSterm
	
	integer :: iCluster, iLayerCluster	
	real*8 :: nParticlesPerCluster, particleCalcite, particleOpal, pressureLimCalc, &
		expFactorCalc, fracDissCalc, dissCalcite, dissOpal, volOpal, volCalc, solidVolDiss, &
		dissolutionRateOpal, expFactorOpal, fracDissOpal
	logical :: carbonZero, tepZero, calciteZero, opalZero, clayZero
	
	! ------ Calculate pressure-based factor for calcite dissolution ------
	pressureLimCalc = CaCO3dissolutionLimFactor(OmegaCalc)
	if (pressureLimCalc > 0d0) then
    	expFactorCalc = EXP(-(dissol_rate_calc/SECONDS_PER_DAY) * pressureLimCalc * timeStep)
    	if (.not. ieee_is_finite(expFactorCalc)) then
    		write(*,*) 'ERROR: expFactorCalc not finite'
    		call WriteStatusAndStop()
		end if
    	fracDissCalc = 1d0 - expFactorCalc
	else
    	fracDissCalc = 0d0
	end if
	fracDissCalc = MAX(0d0, MIN(1d0, fracDissCalc)) ! clamp [0, 1]
	
	! ------ Calculate temperature-dependent factor for opal dissolution ------
	dissolutionRateOpal = dissol_rate_opal_0deg * q10_bSi**((TempC-0d0)/10d0) ! d-1
	expFactorOpal = EXP(-(dissolutionRateOpal/SECONDS_PER_DAY) * timeStep)
	if (.not. ieee_is_finite(fracDissOpal)) then
    	write(*,*) 'ERROR: fracDissOpal not finite'
    	call WriteStatusAndStop()
	end if
	fracDissOpal = 1d0 - expFactorOpal
	fracDissOpal = MAX(0d0, MIN(1d0, fracDissOpal)) ! clamp [0, 1]
		
	do iLayerCluster = 1, nLayerClusters
		iCluster = layerClusterIndices(iLayerCluster)
		if (particle(iCluster)%phase /= 1) cycle ! skip non-particulate clusters	
				
		nParticlesPerCluster = particle(iCluster)%nPxC	 
		if (.not. ieee_is_finite(nParticlesPerCluster) .or. nParticlesPerCluster < 1d0) then
    		write(*,*) 'ERROR: invalid nPxC in mineral dissolution', nParticlesPerCluster
    		call WriteStatusAndStop()
		end if		
		particleCalcite = particle(iCluster)%molesMineral(iCalcite)
		particleOpal = particle(iCluster)%molesMineral(iOpal)
		
		if (particleCalcite <= 0d0 .and. particleOpal <= 0d0) cycle
		
		! ------ Calcite dissolution (abiotic, Omega-limited) ------
		dissCalcite = 0d0			
		if (fracDissCalc > 0d0) then	
			dissCalcite = particleCalcite * fracDissCalc ! mol (in the exponential solution, dissCalcite cannot be > particleCalc)	
			particle(iCluster)%molesMineral(iCalcite) = particleCalcite - dissCalcite
			SMSterm(iDissolCaCO3) = SMSterm(iDissolCaCO3) + dissCalcite*nParticlesPerCluster
			
			! Sanity check: did abiotic dissolution degraded more minerals than existed?
			if (IsExceedingInitialAmount(dissCalcite, particleCalcite, detection_limit_calc)) then
				write(*,*) 'ERROR: more calcite is dissolved than the one available.'
				write(*,*) '  Starting and dissolved material:', particleCalcite, dissCalcite
				write(*,*) '  Depth:', particle(iCluster)%depth
				call WriteStatusAndStop( )
			end if
		end if
					
		! ------ Opal dissolution (temperature-dependent) ------
		dissOpal = 0d0	
		if (particleOpal > 0d0) then
			dissOpal = particleOpal * fracDissOpal ! mol 
			particle(iCluster)%molesMineral(iOpal) = particleOpal - dissOpal
			SMSterm(iDissolOpal) = SMSterm(iDissolOpal) + dissOpal*nParticlesPerCluster
			
			! Sanity check: did abiotic dissolution degraded more minerals than existed?
			if (IsExceedingInitialAmount(dissOpal, particleOpal, detection_limit_opal)) then
				write(*,*) 'ERROR: more opal is dissolved than the one available.'
				write(*,*) '  Starting and dissolved material:', particleOpal, dissOpal
				write(*,*) '  Depth:', particle(iCluster)%depth
				call WriteStatusAndStop( )
			end if
		end if   

		if (dissCalcite + dissOpal <= 0d0) cycle

		nMineralClustersDissolved(iProfile) = nMineralClustersDissolved(iProfile) + 1

		! ------ Shrink ------
		if (shrinkAfterMineralDissolutionScheme == 1 .and. particle(iCluster)%nPpxP >= 2d0 &
			.and. (particle(iCluster)%molesMineral(iOpal) &
			+ particle(iCluster)%molesMineral(iCalcite)) > 0d0) then
			
			volOpal = dissOpal*MOLAR_MASS_OPAL/RHO_OPAL ! cm3
			volCalc = dissCalcite*MOLAR_MASS_CACO3/RHO_CALCITE ! cm3
			solidVolDiss = (volOpal+volCalc)*1d12 ! um3
			if (.not. ieee_is_finite(solidVolDiss) .or. solidVolDiss < 0d0) then
				write(*,*) 'ERROR: invalid dissolved solid volume in mineral dissolution'
    			call WriteStatusAndStop()
			end if
			if (solidVolDiss > 0d0) call ShrinkParticle(particle, nClusters, iCluster, solidVolDiss)
		end if 

		! ------ Recompute particle's derived attributes ------
		call ParticleFractalDimension(particle, nClusters, iCluster, particle(iCluster)%initType)
		call ParticleDryMass(particle, nClusters, iCluster)
		call ParticleMaterialVolume(particle, nClusters, iCluster)
		call ParticleRadius(particle, nClusters, iCluster)
		call ParticlePorosity(particle, nClusters, iCluster)
		call ParticleDensity(particle, nClusters, iCluster, waterRho)
		call ParticleStickiness(particle, nClusters, iCluster)
		call ParticleSettlingVelocity(particle, nClusters, iCluster, waterRho, waterDynVisco)
		
		! ------- Sanity checks ------
		call CheckParticleSanity(particle, nClusters, iCluster, 'mineral dissolution')
		
		! Has the particle become too small? --> flush to dissolution pool (minerals) and solubilisation pool (organics) and mark as empty
		if (2d0*particle(iCluster)%radius < operational_size_poc_min ) then
			SMSterm(iMicrobSolubOrgC) = SMSterm(iMicrobSolubOrgC) + particle(iCluster)%molesOrgC*nParticlesPerCluster
			SMSterm(iMicrobSolubTepC) = SMSterm(iMicrobSolubTepC) + particle(iCluster)%molesTepC*nParticlesPerCluster
			SMSterm(iDissolCaCO3) = SMSterm(iDissolCaCO3) + particle(iCluster)%molesMineral(iCalcite)*nParticlesPerCluster
			SMSterm(iDissolOpal) = SMSterm(iDissolOpal) + particle(iCluster)%molesMineral(iOpal)*nParticlesPerCluster
			SMSterm(iMicrobSolubClay) = SMSterm(iMicrobSolubClay) + particle(iCluster)%molesMineral(iClay)*nParticlesPerCluster
			particle(iCluster)%phase = 3
			cycle ! nothing else to do here
		end if
		
		! Compute total remaining material to derive a tolerance limit
		carbonZero  = IsQuantityEffectivelyZero(particle(iCluster)%molesOrgC, detection_limit_poc)
		tepZero     = IsQuantityEffectivelyZero(particle(iCluster)%molesTepC, detection_limit_poc)
		calciteZero = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iCalcite), detection_limit_calc)
		opalZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iOpal), detection_limit_opal)
		clayZero    = IsQuantityEffectivelyZero(particle(iCluster)%molesMineral(iClay), detection_limit_clay)
		
		! Did the amount of carbon left go beyond the detection limit? --> allocate leftovers to the solubilised pool
		if (calciteZero) then
			SMSterm(iDissolCaCO3) = SMSterm(iDissolCaCO3) + particle(iCluster)%molesMineral(iCalcite)*nParticlesPerCluster
			particle(iCluster)%molesMineral(iCalcite) = 0d0	
		end if
		if (opalZero) then
			SMSterm(iDissolOpal) = SMSterm(iDissolOpal) + particle(iCluster)%molesMineral(iOpal)*nParticlesPerCluster
			particle(iCluster)%molesMineral(iOpal) = 0d0
		end if
	
		! Has the particle been emptied after organics disappeared? --> to empty phase
		if (carbonZero .and. tepZero .and. calciteZero .and. opalZero .and. clayZero) then 
			particle(iCluster)%phase = 3
			cycle
		end if
			
	end do

end subroutine CalciumCarbonateAndOpalDissolution

! ========================================================================================

function CaCO3dissolutionLimFactor(OmegaCalc) 

	!-------------------------------------------------------------------------------------
	! Calculate whether dissolution is possible or not according to Omega. 
	!	If Omega >= 1 --> CO3ion >= CO3sat --> the water is supersaturated with respect to 
	!	mineral CaCO3 --> no dissolution occurs (surface waters).
	!	If Omega < 1 --> the water is undersaturated with respect to mineral CaCO3 -->
	!   dissolution occurs (deep waters)
	!-------------------------------------------------------------------------------------
	
	real*8 :: CaCO3dissolutionLimFactor
	real*8, intent(in) :: OmegaCalc

	if (OmegaCalc >= 1d0) then 
		CaCO3dissolutionLimFactor = 0d0 ! precipitation 
	else
		CaCO3dissolutionLimFactor = (1d0 - OmegaCalc)**order_react_dissol_calc ! dissolution occurs
	end if

end function CaCO3dissolutionLimFactor

! ========================================================================================

end module mineraldissolution