module primaryproduction

! ----------------------------------------------------------------------------------------
! This module handles functions that calculate phytoplankton-related quantities.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use particlestructure, only: lagrangianStateVars
use modelconstants, only: WATT_TO_PHOTON, SECONDS_PER_DAY, MOLAR_MASS_CARBON
use modelparameters, only: nDepthLayers, iEuphoticDepth, maxNumAuxTerms, nPfts, &
	growth_rate_max_0deg_pft, growth_rate_max_phyto, min_par_for_photosynthesis, q10_phyto, &
	k_NO3_pft, k_PO4_pft, k_Si_diat, alpha_chl_spec_pft, Chl2C_ratio, choiceIsPhytoGrowthLim, &
	choicePhytoIrradLimFunc, phyto_life_span_factor, phyto_max_allowed_life_days, iAvgProbDiat, &
	iAvgProbFlagel, iAvgProbCocco, iAvgProbPico, iProfProbDiat, iProfProbFlagel, iProfProbCocco, &
	iProfProbPico
use sanitychecks, only: WriteStatusAndStop
use findfunctions, only: FindDepthLayerIndex
use waterphysicsandlight, only: ParProfile

implicit none
private
public :: EuphoticLayerDepth, PhytoGrowthRate, PhytoSeedingProbability, KillLivingPhytoplanktonCells

contains

! ========================================================================================

subroutine PhytoSeedingProbability(probPftMean, probPftProfile, muPftProfile, validSeedingDepthIdx, &
	nDepthLayersToSeedPhyto, nLocalDepthLayers, auxTerm, auxCount)

	integer, intent(in) :: nDepthLayersToSeedPhyto, nLocalDepthLayers
	real*8, dimension(nPfts,nLocalDepthLayers), intent(in) :: muPftProfile
	logical, dimension(nLocalDepthLayers), intent(in) :: validSeedingDepthIdx
	real*8, dimension(nPfts), intent(out) :: probPftMean
	real*8, dimension(nDepthLayersToSeedPhyto,nPfts), intent(out) :: probPftProfile
	real*8, dimension(maxNumAuxTerms,nDepthLayers), intent(inout) :: auxTerm, auxCount
	
	integer :: iPft, iDl, iValid
	integer, dimension(nPfts), parameter :: pftAvgProbToAux = (/ iAvgProbDiat, iAvgProbFlagel, iAvgProbCocco, iAvgProbPico /)
	integer, dimension(nPfts), parameter :: pftDepthProbToAux = (/ iProfProbDiat, iProfProbFlagel, iProfProbCocco, iProfProbPico /)
	real*8 :: muTotal
	real*8, parameter :: tol_prob = 1d-4
	real*8, dimension(nPfts) :: muPftSum

	! ------ Initialise output ------
	probPftMean(:) = 0d0
    probPftProfile(:,:) = 0d0
    	
	! ------ Sum growth rates across valid depths for each PFT ------
    do iPft = 1, nPfts
        muPftSum(iPft) = SUM(muPftProfile(iPft,:), mask=validSeedingDepthIdx) ! one value per PFT per time step
    end do
    
    ! ------ Total growth rate across all PFTs and depths ------
    muTotal = SUM(muPftSum) 	
	if (.not. ieee_is_finite(muTotal) .or. muTotal < 0d0) then
   		write(*,*) 'ERROR: in PhytoSeedingProbability, muTotal contains NaN or Inf or is negative', muTotal
    	call WriteStatusAndStop( )
	end if
	if (muTotal == 0d0) return

	! ------ Calculate average probability for PFTs across the water column ------
	probPftMean = muPftSum(:)/muTotal ! normalised probability for each PFT
	if (ABS(SUM(probPftMean) - 1d0) > tol_prob) then
    	write(*,*) 'WARNING: probPftMean does not sum to 1:', SUM(probPftMean), probPftMean
	end if

	! ------ Calculate probability for each depth-layer per PFT ------
	do iPft = 1, nPfts
		if (muPftSum(iPft) <= 0d0) cycle
		iValid = 0
        do iDl = 1, nDepthLayersToSeedPhyto
            if (.not. validSeedingDepthIdx(iDl)) cycle
            iValid = iValid + 1
            probPftProfile(iValid,iPft) = muPftProfile(iPft,iDl) / muPftSum(iPft)
        end do
        ! Sanity check: for each PFT separately, the probabilities across depth sum to 1
		if (ABS(SUM(probPftProfile(:,iPft)) - 1d0) > tol_prob) then
   			write(*,*) 'WARNING: probPftProfile does not sum to 1 for PFT', iPft, SUM(probPftProfile(:,iPft))
		end if
	end do
	
	! ------ Store for diagnostics ------ 
	do iPft = 1, nPfts
    	auxTerm(pftAvgProbToAux(iPft),1) = auxTerm(pftAvgProbToAux(iPft),1) + probPftMean(iPft)
        auxCount(pftAvgProbToAux(iPft),1) = auxCount(pftAvgProbToAux(iPft),1) + 1d0
    end do	
	do iDl = 1, nDepthLayersToSeedPhyto
    	do iPft = 1, nPfts
        	auxTerm(pftDepthProbToAux(iPft),iDl) = auxTerm(pftDepthProbToAux(iPft),iDl) + probPftProfile(iDl,iPft)
        	auxCount(pftDepthProbToAux(iPft),iDl) = auxCount(pftDepthProbToAux(iPft),iDl) + 1d0
    	end do
	end do

end subroutine PhytoSeedingProbability

! ========================================================================================

subroutine PhytoGrowthRate(muPftProfile, validSeedingDepthIdx, nDepthLayersToSeedPhyto, &
	nLocalDepthLayers, parz, NO3, PO4, SiOH4, TempC)

	integer, intent(in) :: nDepthLayersToSeedPhyto, nLocalDepthLayers
	real*8, dimension(nDepthLayersToSeedPhyto), intent(in) :: NO3, PO4, SiOH4, TempC, parz
	logical, dimension(nLocalDepthLayers), intent(in) :: validSeedingDepthIdx
	real*8, dimension(nPfts,nLocalDepthLayers), intent(inout) :: muPftProfile

	integer :: iDl, iPft
	real*8 :: tempFunc, muMax0deg, muT, nutLim, chlToCratio, initSlopePIcurve, irrLim, muLimFactor
	real*8, dimension(nDepthLayersToSeedPhyto) :: parzumol
	
	parzumol(:) = parz(:)*WATT_TO_PHOTON ! W m-2 --> umol photons m-2 s-1

	do iDl = 1, nDepthLayersToSeedPhyto
		if (.not. validSeedingDepthIdx(iDl)) cycle   ! skip non-euphotic depths
        if (NO3(iDl) <= 0d0) cycle                   ! skip if no nitrate

		do iPft = 1, nPfts
			! Nutrient-and-temperature limiting growth rate
			if (.not. ieee_is_finite(TempC(iDl))) cycle
			tempFunc = PhytoTempFunc(iPft, TempC(iDl))
			muMax0deg = growth_rate_max_0deg_pft(iPft)
			muT = muMax0deg * tempFunc ! d-1

			select case (choiceIsPhytoGrowthLim)
			case (1) ! there is growth limitation
				nutLim = PhytoNutLim(iPft, NO3(iDl), PO4(iDl), SiOH4(iDl))
				chlToCratio = PhytoChlToCratio(iPft, muT, parzumol(iDl)) ! g Chl (g C)-1
				initSlopePIcurve = chlToCratio*alpha_chl_spec_pft(iPft) ! m2 (umol photons)-1
				if (initSlopePIcurve <= 0d0) cycle
				irrLim = PhytoIrrLim(muT, initSlopePIcurve, parzumol(iDl))

			case (0) ! there is no growth limitation
				nutLim = 1d0
				irrLim = 1d0
			end select

			! Growth rate
			muLimFactor = nutLim*irrLim
			muPftProfile(iPft,iDl) = MIN((muT*muLimFactor), growth_rate_max_phyto) ! divisions d-1
		end do
	end do
	
	! Sanity check 
	if (.not. all(ieee_is_finite(muPftProfile)) .or. any(muPftProfile < 0d0)) then
   		write(*,*) 'ERROR: in PhytoGrowthRate, muPftProfile contains NaN or Inf or negative values'
    	call WriteStatusAndStop( )
	end if
	
end subroutine PhytoGrowthRate
	    
! ========================================================================================

subroutine EuphoticLayerDepth(zeu, parz, nDepthLayersToSeedPhyto, validSeedingDepthIdx, auxTerm, &
	auxCount, Chla, PAR0, MLD, nLocalDepthLayers, zmid)

	integer, intent(in) :: nLocalDepthLayers
	real*8, intent(in) :: Chla, PAR0, MLD	
	real*8, dimension(nLocalDepthLayers), intent(in) :: zmid
	real*8, dimension(maxNumAuxTerms), intent(inout) :: auxTerm, auxCount
	integer, intent(out) :: nDepthLayersToSeedPhyto
	real*8, intent(out) :: zeu
	real*8, dimension(nLocalDepthLayers), intent(out) :: parz
	logical, dimension(nLocalDepthLayers), intent(out) :: validSeedingDepthIdx
	
	integer :: iDl, lastTrue
	logical, dimension(nLocalDepthLayers) :: tmpValidSeedingDepthIdx
	
	! Compute PAR at each depth
	call ParProfile(parz(:), nLocalDepthLayers, zmid(:), PAR0, Chla, MLD) ! --in waterphysicsandoptics.F90

	! Identify valid depths for phytoplankton seeding
	if (PAR0 <= 0d0) then
    	tmpValidSeedingDepthIdx(:) = .false.
	else
    	tmpValidSeedingDepthIdx = (parz(:)/PAR0) >= min_par_for_photosynthesis
	end if

	! Find deepest (last) index where condition holds
	lastTrue = 0
	do iDl = nLocalDepthLayers, 1, -1
		if (tmpValidSeedingDepthIdx(iDl)) then
			lastTrue = iDl
			exit
		end if
	end do
	
	! Build the final valid mask: contiguous from surface (1) down to lastTrue
	validSeedingDepthIdx(:) = .false.
	if (lastTrue > 0) validSeedingDepthIdx(1:lastTrue) = .true.
	if (lastTrue > 0) then
  		nDepthLayersToSeedPhyto = lastTrue
  		zeu = zmid(lastTrue)
	else
  		nDepthLayersToSeedPhyto = 0
  		zeu = 0d0
	end if
	
	! Sanity checks
	if (.not. ieee_is_finite(zeu)) then
  		write(*,*) 'ERROR: invalid zeu:', zeu
  		call WriteStatusAndStop( )
	end if
	if (.not. all(ieee_is_finite(parz))) then
    	write(*,*) 'ERROR: parz contains NaN or Inf', parz
    	call WriteStatusAndStop( )
	end if
	
	! Store zeu value for diagnostics 
	if (nDepthLayersToSeedPhyto > 0) then	
		auxTerm(iEuphoticDepth) = auxTerm(iEuphoticDepth) + zeu
		auxCount(iEuphoticDepth) = auxCount(iEuphoticDepth) + 1d0
	else
		auxTerm(iEuphoticDepth) = auxTerm(iEuphoticDepth) + 0d0
		auxCount(iEuphoticDepth) = auxCount(iEuphoticDepth) + 1d0
	end if

end subroutine EuphoticLayerDepth

! ========================================================================================

function PhytoTempFunc(iPft, TempC)

	real*8 :: PhytoTempFunc
	integer, intent(in) :: iPft
	real*8, intent(in) :: TempC
	
	if (.not. ieee_is_finite(TempC)) then
    	write(*,*) 'ERROR: invalid TempC in PhytoTempFunc'
    	call WriteStatusAndStop()
	end if

	! Eppley function (Eppley 1972)
	PhytoTempFunc = q10_phyto(iPft)**((TempC-0d0)/10d0) ! =exp(K_EPPLEY*TempC)

end function PhytoTempFunc

! ========================================================================================

function PhytoNutLim(iPft, NO3, PO4, SiOH4)

	real*8 :: PhytoNutLim
	integer, intent(in) :: iPft
	real*8, intent(in) :: NO3, PO4, SiOH4
	real*8 :: fNO3, fPO4, fSi
	
	fNO3 = MAX(NO3, 0d0) / (MAX(NO3, 0d0) + k_NO3_pft(iPft))
    fPO4 = MAX(PO4, 0d0) / (MAX(PO4, 0d0) + k_PO4_pft(iPft))

	! Apply Michaelis-Menten kinetics based on the PFT and Liebig's law of the minimum
    if (iPft == 1) then
        ! Diatoms require NO3, PO4, and SiOH4
        fSi = MAX(SiOH4, 0d0) / (MAX(SiOH4, 0d0) + k_Si_diat)
        PhytoNutLim = MIN(fNO3, fPO4, fSi)
    else
        ! Other PFTs require only NO3 and PO4
        PhytoNutLim = MIN(fNO3, fPO4)
    end if

end function PhytoNutLim

! ========================================================================================

function PhytoIrrLim(muT, initSlopePIcurve, PAR)

	real*8 :: PhytoIrrLim
	real*8, intent(in) :: muT, initSlopePIcurve, PAR
	real*8 :: kPAR, photoinhibitionCoeff, PARinh

	! Sanity checks
	if (.not. ieee_is_finite(muT) .or. muT <= 0d0) then
		PhytoIrrLim = 0d0
		return
	end if
	if (.not. ieee_is_finite(initSlopePIcurve) .or. initSlopePIcurve <= 0d0) then
		PhytoIrrLim = 0d0
		return
	end if
	if (.not. ieee_is_finite(PAR) .or. PAR <= 0d0) then
		PhytoIrrLim = 0d0
		return
	end if

	kPAR = (muT/SECONDS_PER_DAY) / initSlopePIcurve ! umol photons m-2 s-1 

	! Choice of an empirical function that establishes a dependency of photosynthesis with 
	! light (i.e., PI curve). These are analytical expressions
	select case (choicePhytoIrradLimFunc)
	case (1) ! Michaelis-Menten-Monod hyperbolic saturation model, Baly (1935)
		PhytoIrrLim = PAR / (kPAR+PAR)

	case (2) ! Exponential saturation model, Webb et al. (1974)
		PhytoIrrLim = 1d0 - EXP(-PAR/kPAR)

	case (3) ! Hyperbolic saturation model (Smith function), Platt & Jassby (1976)
		PhytoIrrLim = PAR / (sqrt(kPAR**2 + PAR**2)) ! equivalent to irrLim = V/muT, where V = alpha*PARz*muT/ (sqrt(muT^2 + (alpha*PARz)^2))

	case (4) ! Exponential function with photoinhibition at hight PAR, Steele (1962)
		PhytoIrrLim = (PAR/kPAR) * EXP(1d0 - (PAR/kPAR))

	case (5) ! Hyperbolic tangent saturation model, Platt & Jassby (1976)
		PhytoIrrLim = TANH(PAR/kPAR)

	case (6) ! Exponential saturation and photoinhibition model, Platt et al. (1980)
		photoinhibitionCoeff = 5.963d-15 ! umol photons m-2 s-1 (mem model)
		PARinh = (muT/SECONDS_PER_DAY) / photoinhibitionCoeff
		PhytoIrrLim = (1d0 - EXP(-PAR/kPAR)) * EXP(-PAR/PARinh)

	end select
	
	! In the model MEDUSA, the equivalent expressions are:	
	! 	kPAR = muT / initSlopePIcurve (Ik = muT / alpha)
	! 	PhytoIrrLim = V/muT, V is the expression J in MEDUSA

end function PhytoIrrLim

! ========================================================================================

function PhytoChlToCratio(iPft, muT, PAR)

	real*8 :: PhytoChlToCratio
	integer, intent(in) :: iPft
	real*8, intent(in) :: muT, PAR

	if (.not. ieee_is_finite(muT) .or. muT <= 0d0) then
    	PhytoChlToCratio = Chl2C_ratio
    	return
	end if

	! From Geider et al. (1997)	
	PhytoChlToCratio = Chl2C_ratio / &
		(1d0 + (Chl2C_ratio*alpha_chl_spec_pft(iPft)*PAR / (2d0*muT/SECONDS_PER_DAY))) ! g Chl (g C)-1

end function PhytoChlToCratio

! ========================================================================================

subroutine KillLivingPhytoplanktonCells(particle, nClusters, iCluster, zeu, muPftProfile, &
	nLocalDepthLayers, zbot, ztop, iTimeStep, timeStepDays)
	
	!-------------------------------------------------------------------------------------
	! Kill living cells if they sink below zeu or exceed maximum age.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nClusters, nLocalDepthLayers, iCluster, iTimeStep
	real*8, intent(in) :: zeu, timeStepDays
	real*8, dimension(nPfts,nLocalDepthLayers), intent(in) :: muPftProfile
	real*8, dimension(nLocalDepthLayers), intent(in) :: zbot, ztop
	type(lagrangianStateVars), dimension(nClusters), intent(inout) :: particle
	
	integer :: iDepthLayer, pftType
	real*8 :: cellAgeDays, muPftLocal, doublingTime, maxLifeDays
	
	if (particle(iCluster)%living /= 1) return
	
	! PFT
	pftType = particle(iCluster)%initPft
		
	! Calculate cell age
	cellAgeDays = dble(iTimeStep - particle(iCluster)%tstepCreat)*timeStepDays
	if (.not. ieee_is_finite(cellAgeDays)) then
		write(*,*) 'ERROR in KillLivingPhytoplanktonCells: invalid cellAgeDays =', cellAgeDays, ' for pft=', pftType
		call WriteStatusAndStop( )
	end if
	
	! Extract depth layer
	iDepthLayer = FindDepthLayerIndex(particle(iCluster)%depth, nLocalDepthLayers, ztop(:), zbot(:))
	if (iDepthLayer < 1 .or. iDepthLayer > nLocalDepthLayers) then	
		write(*,*) 'ERROR in KillLivingPhytoplanktonCells: particle depth not in any layer: depth=', particle(iCluster)%depth
		call WriteStatusAndStop( )
	end if
	
	! Local growth rate
	muPftLocal = muPftProfile(pftType,iDepthLayer)
	if (.not. ieee_is_finite(muPftLocal) .or. muPftLocal < 0d0) then
		write(*,*) 'ERROR in KillLivingPhytoplanktonCells: invalid local growth rate muLocal=', muPftLocal, ' for pft=', pftType, ' layer=', iDepthLayer
		call WriteStatusAndStop( )
	end if
	
	! Decide life span
	if (muPftLocal > 0d0) then
		doublingTime = LOG(2d0)/muPftLocal ! days
		maxLifeDays = phyto_life_span_factor * doublingTime
		if (maxLifeDays > phyto_max_allowed_life_days) maxLifeDays = phyto_max_allowed_life_days ! clamp life to avoid ridiculously large lives for tiny mu
	else	
		! no growth potential -> make a conservative short life (starvation)
		maxLifeDays = 1d0 ! tune: 1 day for fully starving cells
	end if
	
	! Mortality rule: death by age OR by being deeper than zeu
	if (cellAgeDays >= maxLifeDays .or. particle(iCluster)%depth > zeu) then
		particle(iCluster)%living = 0
		particle(iCluster)%initType = 4
	end if

end subroutine KillLivingPhytoplanktonCells

! ========================================================================================

end module primaryproduction