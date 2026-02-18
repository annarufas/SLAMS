#include "blockdefinitions.h"

module initialisation

! ----------------------------------------------------------------------------------------
! This module is concerned with initialising the model: creating and initialising the 
! particle array, creating and initialising the model counters, reading the namelist 
! parameters, creating the physical grid, loading the environmental (forcing) data and
! creating and initialising the model output collection arrays (SMS and auxiliary terms, 
! sediment traps and imaging systems)
! ----------------------------------------------------------------------------------------

use modelcounters, only: InitialiseCounters
use modelparameters, only: GetDepthsForParticleDataCollection, EstablishParticleSizeAndVeloCategories, &
	InitialiseRuntimeParameters, EstablishMesozooSizeCategoriesAndSpectraVariables
use modelgrid, only: AreaGrid, DepthLayerGrid
use modelforcingdata, only: LoadBinForcingData, IniForcingData
use modeleulerianvariables, only: EulerianVariables
use modelparticlecollection, only: ParticleCollectionVariables
use sanitychecks, only: WriteStatus

implicit none
public :: InitialiseModel

contains

! ========================================================================================

subroutine InitialiseModel(nProfiles)

	!-------------------------------------------------------------------------------------
	! This is for global, rather than local, configurations.
	!-------------------------------------------------------------------------------------

	integer, intent(in) :: nProfiles

	call WriteStatus('in progress')
	call InitialiseCounters(nProfiles) 							! --in modelcounters.F90
	call GetDepthsForParticleDataCollection( ) 					! --in modelparameters.F90
	call EstablishParticleSizeAndVeloCategories( ) 				! --in modelparameters.F90
	write(*,*) 'About to read namelist'
	call InitialiseRuntimeParameters( ) 						! --in modelparameters.F90
	write(*,*) 'Namelist read OK'
    call AreaGrid(nProfiles)                                    ! --in modelgrid.F90
#ifndef BLOCK_TMM_SLAMS
	call DepthLayerGrid( ) 										! --in modelgrid.F90
	call LoadBinForcingData( ) 									! --in modelforcingdata.F90
#endif
    call IniForcingData(nProfiles) 								! --in modelforcingdata.F90
	call EulerianVariables(nProfiles) 							! --in modeleulerianvariables.F90
	call ParticleCollectionVariables(nProfiles) 				! --in modelparticlecollection.F90
	call EstablishMesoZooSizeCategoriesAndSpectraVariables( ) 	! --in modelparameters.F90

end subroutine InitialiseModel

! ========================================================================================

end module initialisation