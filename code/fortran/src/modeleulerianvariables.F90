module modeleulerianvariables

! ----------------------------------------------------------------------------------------
! This module is concerned with carrying the arrays that save the model SMS terms. 
! ----------------------------------------------------------------------------------------

use modelparameters, only: nDepthLayers, maxNumSmsTerms, maxNumTracers, maxNumAuxTerms, &
	iPrimProdOrgC, iPrimProdCaCO3, iPrimProdOpal, iMicrobSolubOrgC, iMicrobSolubTepC, &
	iMicrobSolubCaCO3, iMicrobSolubOpal, iMicrobSolubClay, iZooSolubOrgC, iZooSolubTepC, &
	iZooSolubCaCO3, iZooSolubOpal, iZooSolubClay, iPhotoTepC, iZooRespOrgC, iZooRespTepC, &
	iZooDissolCaCO3, iZooExcretOrgC, iZooExcretTepC, iMicrobRespOrgC, iMicrobRespTepC, &
	iDissolCaCO3, iDissolOpal, C2P_ratio, C2N_ratio	

implicit none

real*8, dimension(:,:,:), allocatable :: SMSterm
real*8, dimension(:,:,:), allocatable :: accumSMS
real*8, dimension(:,:,:), allocatable, target :: tracerSMSterms
real*8, dimension(:,:,:), allocatable :: avgSMSflux
real*8, dimension(:,:), allocatable :: avgSMSfluxIntegrated

real*8, dimension(:,:,:), allocatable :: auxTerm
real*8, dimension(:,:,:), allocatable :: auxCount

contains

! ========================================================================================

subroutine EulerianVariables( nProfiles )

	integer, intent(in) :: nProfiles

	allocate(SMSterm(maxNumSmsTerms,nDepthLayers,nProfiles))
	allocate(accumSMS(maxNumSmsTerms,nDepthLayers,nProfiles))
	allocate(tracerSMSterms(nDepthLayers,maxNumTracers,nProfiles))
	allocate(avgSMSflux(maxNumSmsTerms,nDepthLayers,nProfiles))
	allocate(avgSMSfluxIntegrated(maxNumSmsTerms,nProfiles))

	allocate(auxTerm(maxNumAuxTerms,nDepthLayers,nProfiles))
	allocate(auxCount(maxNumAuxTerms,nDepthLayers,nProfiles))
	
	SMSterm(:,:,:)            = 0d0
	accumSMS(:,:,:)           = 0d0
	tracerSMSterms(:,:,:)     = 0d0
	avgSMSflux(:,:,:)         = 0d0
	avgSMSfluxIntegrated(:,:) = 0d0
	
	auxTerm(:,:,:)            = 0d0
	auxCount(:,:,:)           = 0d0

end subroutine EulerianVariables

! ========================================================================================

subroutine SourcesMinusSinksTermInformation( nLocalDepthLayers, iProfile )

	integer, intent(in) :: nLocalDepthLayers, iProfile
    
    integer :: iz

!	This function is called every time step.
    
	accumSMS(:,1:nLocalDepthLayers,iProfile) = accumSMS(:,1:nLocalDepthLayers,iProfile) &
		+ SMSterm(:,1:nLocalDepthLayers,iProfile)

    do iz = 1, nLocalDepthLayers
    
!  	  DOC, mol
	  tracerSMSterms(iz,1,iProfile) = &
		  + SMSterm(iMicrobSolubOrgC,iz,iProfile) &
		  + SMSterm(iMicrobSolubTepC,iz,iProfile) &
		  + SMSterm(iZooSolubOrgC,iz,iProfile) &
		  + SMSterm(iZooSolubTepC,iz,iProfile) &
		  + SMSterm(iZooExcretOrgC,iz,iProfile) &
		  + SMSterm(iZooExcretTepC,iz,iProfile) &
		  + SMSterm(iPhotoTepC,iz,iProfile)
		  		  
!  	  DIC (CO2 + CO3=), mol
	  tracerSMSterms(iz,2,iProfile) = &
		  - SMSterm(iPrimProdOrgC,iz,iProfile) &
		  - SMSterm(iPrimProdCaCO3,iz,iProfile) &
		  + SMSterm(iMicrobRespOrgC,iz,iProfile) &
		  + SMSterm(iMicrobRespTepC,iz,iProfile) &
		  + SMSterm(iZooRespOrgC,iz,iProfile) &
		  + SMSterm(iZooRespTepC,iz,iProfile) &
		  + SMSterm(iDissolCaCO3,iz,iProfile) &
		  + SMSterm(iZooDissolCaCO3,iz,iProfile) &
		  + SMSterm(iMicrobSolubCaCO3,iz,iProfile) &
		  + SMSterm(iZooSolubCaCO3,iz,iProfile)
		
!	  Silicic acid, mol		
	  tracerSMSterms(iz,3,iProfile) = &
		  - SMSterm(iPrimProdOpal,iz,iProfile) &
		  + SMSterm(iDissolOpal,iz,iProfile) &
		  + SMSterm(iMicrobSolubOpal,iz,iProfile) &
		  + SMSterm(iZooSolubOpal,iz,iProfile)
		  		
!     Phosphate, mol
	  tracerSMSterms(iz,4,iProfile) = &
		  tracerSMSterms(iz,1,iProfile)/C2P_ratio
	
! 	  Nitrate, mol
	  tracerSMSterms(iz,5,iProfile) = &
		  tracerSMSterms(iz,1,iProfile)/C2N_ratio 

! 	  ALK = [HCO3-] + 2[CO32-] + [OH-] - [H+] + 2[PO43-]
	  tracerSMSterms(iz,6,iProfile) = &
		  - 2d0*SMSterm(iPrimProdCaCO3,iz,iProfile) &
		  - 2d0*SMSterm(iPrimProdOrgC,iz,iProfile)/C2N_ratio &
		  + 2d0*SMSterm(iDissolCaCO3,iz,iProfile) &
		  + 2d0*SMSterm(iZooDissolCaCO3,iz,iProfile) &
		  + 2d0*SMSterm(iMicrobSolubCaCO3,iz,iProfile) &
		  + 2d0*SMSterm(iZooSolubCaCO3,iz,iProfile)
		  
    end do 				
 	  			  
end subroutine SourcesMinusSinksTermInformation

! ========================================================================================

end module modeleulerianvariables