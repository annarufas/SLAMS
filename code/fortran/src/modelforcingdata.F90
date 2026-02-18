module modelforcingdata

! ----------------------------------------------------------------------------------------
! This module handles the reading of model forcing data from input files.
! ----------------------------------------------------------------------------------------

use modelparameters, only: nDepthLayers, filenameNPP, filenameChla, filenameMLD, &
	filenamePAR0, filenameAeolClay, filenameNO3, filenamePO4, filenameSiOH4, filenameMesoZoo, &
	filenameOmegaCalc, filenameTempC, filenameO2, filenameRho, filenameDynVisco, nTimeStepsYear

implicit none
external :: read_r8_field
public

save ! ensure module variables persist between calls
real*8, dimension(:), allocatable   :: NPP, Chla, MLD, PAR0, AeolClay
real*8, dimension(:,:), allocatable :: NO3, PO4, SiOH4, MesoZoo, OmegaCalc, TempC, O2, Rho, DynVisco
real*8, dimension(:), allocatable   :: accumAeolClay

contains

! ========================================================================================

subroutine LoadBinForcingData( )

	allocate(NPP(nTimeStepsYear)) 					! Net primary production, mol C m-2 s-1
	allocate(Chla(nTimeStepsYear)) 					! Chlorophyll a concentration, mg m-3
	allocate(MLD(nTimeStepsYear))					! Mixed layer depth, m
	allocate(PAR0(nTimeStepsYear)) 					! Photosynthetic active radiation at the surface ocean, W m-2 (= J s-1 m-2)
	allocate(AeolClay(nTimeStepsYear)) 				! Aeolian clay deposition, g clay m-2 s-1

	allocate(NO3(nDepthLayers,nTimeStepsYear)) 		! Nitrate concentration, mmol m-3
	allocate(PO4(nDepthLayers,nTimeStepsYear)) 		! Phosphate concentration, mmol m-3
	allocate(SiOH4(nDepthLayers,nTimeStepsYear)) 	! Silicic acid concentration, mmol m-3
	allocate(MesoZoo(nDepthLayers,nTimeStepsYear)) 	! Mesozooplankton concentration, mg C m-3
	allocate(OmegaCalc(nDepthLayers,nTimeStepsYear))! Omega calcite
	allocate(TempC(nDepthLayers,nTimeStepsYear)) 	! Temperature, ºC
	allocate(O2(nDepthLayers,nTimeStepsYear))		! Oxygen concentration, mL L-1 
	allocate(Rho(nDepthLayers,nTimeStepsYear))	 	! Seawater density, g cm-3
	allocate(DynVisco(nDepthLayers,nTimeStepsYear))	! Dynamic viscosity of seawater, g cm-1 s-1							
                                         
	call read_r8_field(nTimeStepsYear,1,1,NPP,filenameNPP)
	call read_r8_field(nTimeStepsYear,1,1,Chla,filenameChla)
	call read_r8_field(nTimeStepsYear,1,1,MLD,filenameMLD) 
	call read_r8_field(nTimeStepsYear,1,1,PAR0,filenamePAR0) 
	call read_r8_field(nTimeStepsYear,1,1,AeolClay,filenameAeolClay) 
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,NO3,filenameNO3)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,PO4,filenamePO4)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,SiOH4,filenameSiOH4)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,MesoZoo,filenameMesoZoo)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,OmegaCalc,filenameOmegaCalc)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,TempC,filenameTempC)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,O2,filenameO2)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,Rho,filenameRho)
	call read_r8_field(nDepthLayers,nTimeStepsYear,1,DynVisco,filenameDynVisco)
			
end subroutine LoadBinForcingData

! ========================================================================================

subroutine IniForcingData( nProfiles )

! Initialise additional arrays required for forcing data computations

	integer, intent(in) :: nProfiles

	allocate(accumAeolClay(nProfiles))
	
	accumAeolClay(:) = 0d0

end subroutine IniForcingData

! ========================================================================================

end module modelforcingdata