! ========================================================================================
! 										SLAMS-2.0
! 	a Stochastic, Lagrangian Aggregate Model of Sinking particles in the ocean, v. 2.0
! ========================================================================================

#include "blockdefinitions.h"

program BiologicalCarbonPump

! ----------------- Declaration of modules used in main -----------------
use particlestructure
use modelcounters, only: iTimeStepYear, iYear, nZooDeadClustersPerProfile
use modelparameters, only: nYears, nDepthLayers, maxNumClusters, maxNumClustersPerProfile, &
	nTimeStepsPassedToShiftClusters, maxNumTimeSteps, nTimeStepsYear, nTimeStepsDay, &
	maxNumTracers, nPfts	
use initialisation
use modelgrid, only: gridCellArea, ztop, zbot, zmid
use modelforcingdata, only: NPP, Chla, PAR0, MLD, AeolClay, NO3, PO4, SiOH4, MesoZoo, DynVisco, &
	OmegaCalc, TempC, O2, Rho
use modeleulerianvariables, only: SMSterm, tracerSMSterms, SourcesMinusSinksTermInformation, &
	auxTerm, auxCount 
use modelparticlecollection, only: nSedTrapDeployDepths, nImagingDeployDepths, &
	sedTrapDeployDepths, sedTrap, sedTrapSf, attsInSizeClass, attsInVeloClass, attsInVolumeClass, &
	attsInSizeClassSf, attsInVeloClassSf, attsInMainType, attsInMainTypeSf
use sanitychecks, only: WriteStatus, WriteStatusAndStop
use findfunctions, only: FindActiveClusters, FindClustersInDepthLayer, FindPhase1clusters, &
	FindDepthLayerMidpointDepth
use calcparticleattributes, only: CompactParticleArray
use waterphysicsandlight, only: TurbulentKineticEnergyDissipationRateDepthProfile
use heterotrophicmetabolism, only: ComputeBianchiDvmBounds, ZooplanktonDielVerticalMigrationDepth
use modeloutput, only: GetParticleAverageAttributes, GetParticleAverageAttributesByVolumeClass, &
	WriteSnapshotsOfParticleAttributes, WriteInstantaneousSnapshots, WriteModelClosureInformation
use wrappers, only: ProductionOfSurfaceParticles, ParticleProcessing, GravitationalSettling, &
	WritePeriodicOutput
use timer, only: StepTimer, IniStepTimer, UpdateStepTimer
use theseed, only: inititialiseRandomSeed90, inititialiseRandomSeed95

implicit none

! ----------------- Declaration of variables used in main -----------------
integer :: iTimeStep, iForcingStep, iDepthLayer, iCluster, iActiveCluster, iProfile, &
	iLastLocus, nProfiles, nClusters, nLocalDepthLayers, nSurfLayersBianchi, nMesoLayersBianchi, &
	nLayerClusters, nActiveClusters, nPhase1clusters, iLastClusterPre
integer, dimension(:), allocatable :: layerClusterIndices, activeClusterIndices, phase1clusterIndices
integer, dimension(4) :: dvmBoundsBianchi
real*8 :: startTimer, stopTimer, ztopLayer, zbotLayer, zmidLayer, gridCellVolume, localSeafloorDepth, &
	tkeLayer, dvmDepthNightUpper, dvmDepthNightLower, dvmDepthDayUpper, dvmDepthDayLower, zeu
real*8, dimension(:), allocatable :: tkeDepthProfile
real*8, dimension(:,:), allocatable :: muPftProfile
type(lagrangianStateVars), dimension(:), allocatable :: particle

integer, parameter :: nMonths = 12
integer, dimension(nMonths), parameter :: daysInMonth = &
    (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
integer, dimension(nMonths) :: midMonthTimeStep
integer :: m, dayStart, dayMid
logical :: isMidMonth

! ----------------- Declaration of diagnostics flag -----------------
#ifndef BLOCK_TMM_SLAMS
#ifdef BLOCK_DIAGNOSTICS
	type(StepTimer) :: diagTimer
	logical :: doDiagnosticsAvg = .FALSE. ! flag to indicate whether the end of a calendar month has been reached
	integer :: diagsAccumTimeSteps
#endif
	logical :: doDiagnosticsAccum = .FALSE. ! flag to indicate whether diagnostics should be accumulated
#endif

! ----------------------------------------------------------------------------------------
! ----------------------------------MODEL INITIALISATION----------------------------------
! ----------------------------------------------------------------------------------------

! ------ Initialise random seed for unique sequences ------
#ifdef BLOCK_NEW_RANDOM_SEED
	call inititialiseRandomSeed95()
#endif

! ------ Run this for a timer ------
dayStart = 1
do m = 1, nMonths
    dayMid = dayStart + daysInMonth(m)/2   ! integer mid-day
    midMonthTimeStep(m) = (dayMid - 1) * nTimeStepsDay + 1
    ! "+1" picks the first timestep of the day (since 3 per day)
    dayStart = dayStart + daysInMonth(m)
end do
write(*,*) 'midMonthTimeStep:', midMonthTimeStep

! ------ Initialise model framework ------
nProfiles = 1 ! all the processed performed by this program happen in one water column profile
call InitialiseModel(nProfiles) ! --in initialisation.F90

! ------ Set local depth layers and seafloor depth ------
nLocalDepthLayers = nDepthLayers
localSeafloorDepth = zbot(nLocalDepthLayers) ! m
call ComputeBianchiDvmBounds(dvmBoundsBianchi(:), nSurfLayersBianchi, nMesoLayersBianchi, &
	nLocalDepthLayers, zmid(:), localSeafloorDepth) ! --in heterotrophicmetabolism.F90
allocate(muPftProfile(nPfts,nLocalDepthLayers))
allocate(tkeDepthProfile(nLocalDepthLayers))

#ifdef BLOCK_FEW_TIME_STEPS
	maxNumTimeSteps = 21 ! Optional: limit the number of time steps to 21 for testing/debugging
#endif

! ------ Initialise the particle array (Lagrangian variables) ------
! All clusters start with zeroed attributes (ID = 0), meaning they are unused. Active 
! clusters have phase = 1 (water column), while phase = 2 (seafloor) and phase = 3 (emptied) 
! are inactive and await reuse.
nClusters = maxNumClusters
allocate(particle(nClusters))
call InitialiseParticleArray(particle, nClusters) ! --in particlestructure.F90

! ------ Only proceed if the water column depth exceeds 25 m ------
if (localSeafloorDepth <= 25d0) then
	call WriteStatusAndStop('water too shallow')
end if

iLastLocus = 0              ! initialise index for particle array
nActiveClusters = 0         ! initialise nActiveClusters
iProfile = 1                ! set current profile index (used by TMM interface)
call cpu_time( startTimer ) ! start CPU timer

write(*,*)
write(*,*) '=========================== STARTING SLAMS2.0 ==========================='
write(*,*)
write(*,'(A,F8.1,A)') 'The water column is ', localSeafloorDepth, ' m depth'
write(*,'(A,I0,A)') 'The number of depth layers is ', nLocalDepthLayers
write(*,*)

#ifdef BLOCK_DIAGNOSTICS
#ifndef BLOCK_TMM_SLAMS
	call IniStepTimer(0, diagTimer) ! initialise the diagnostic timer
#endif
#endif

! ----------------------------------------------------------------------------------------
! ------------------------------------TIME STEP BEGINS------------------------------------
! ----------------------------------------------------------------------------------------

do iTimeStep = 1, maxNumTimeSteps

	! Update time step counters
	if (iProfile == 1) then
		iTimeStepYear = iTimeStepYear + 1 
	
		! If a year has passed, reset iTimeStepYear and update the year
		if (MOD((iTimeStep-1),nTimeStepsYear) == 0) then
			iTimeStepYear = 1
			iYear = iYear + 1
			write(*,'(A,I0,A)') '************************** Year ', iYear, ' **************************'
		end if
		
		! Get mid-month day
		isMidMonth = .false.
		do m = 1, nMonths
			if (iTimeStepYear == midMonthTimeStep(m)) then
				isMidMonth = .true.
				exit
			end if
		end do			
	end if 

	iForcingStep = iTimeStepYear ! set forcing time step to current iTimeStepYear

#ifdef BLOCK_DIAGNOSTICS
#ifndef BLOCK_TMM_SLAMS
	if (iTimeStep>=diagTimer%startTimeStep) doDiagnosticsAccum = .TRUE.
#endif
#endif

	! Reset TKE and SMS terms for the current profile at the beginning of the time step
	SMSterm(:,:,iProfile) = 0d0
	tracerSMSterms(:,:,iProfile) = 0d0
	tkeDepthProfile(:) = 0d0

	! The TKE surface value changes at every time step and TKE has to be calculated for the 
	! whole profile with that surface value						
	call TurbulentKineticEnergyDissipationRateDepthProfile( tkeDepthProfile(:), &
		nLocalDepthLayers, localSeafloorDepth ) ! --in waterphysicsandoptics.F90

#ifdef BLOCK_TMM_SLAMS
	nExistinglusters = nActiveClusters - maxNumClustersPerProfile
	iLastLocusPre = iLastLocus
#endif

! ----------------------------------------------------------------------------------------
! -------------------------------ADD PARTICLES TO THE ARRAY-------------------------------
! ----------------------------------------------------------------------------------------

	! Seed phytoplankton, TEP and clay.
	call ProductionOfSurfaceParticles(particle, nClusters, iLastLocus, zeu, nLocalDepthLayers, &
		ztop(:), zbot(:), zmid(:), gridCellArea(iProfile), SMSterm(:,:,iProfile), auxTerm(:,:,iProfile), &
		auxCount(:,:,iProfile), muPftProfile(:,:), NPP(iForcingStep), Chla(iForcingStep), &
		PAR0(iForcingStep), MLD(iForcingStep), AeolClay(iForcingStep), NO3(:,iForcingStep), &
		PO4(:,iForcingStep), SiOH4(:,iForcingStep), OmegaCalc(:,iForcingStep), TempC(:,iForcingStep), &
		Rho(:,iForcingStep), DynVisco(:,iForcingStep), iTimeStep, iProfile) ! --in wrappers.F90
		 
#ifdef BLOCK_TMM_SLAMS
	nNewClusters = iLastLocus - iLastLocusPre
#endif

#ifndef BLOCK_TMM_SLAMS
	if (iLastLocus > 0) then
#else
	if (nExistinglusters + nNewClusters > 0) then
#endif
	
#ifdef BLOCK_DIAGNOSTICS		
#ifndef BLOCK_FEW_TIME_STEPS	
	if (doDiagnosticsAccum) then
#endif
	! Image current particles and classify them according to their size and sinking speed.
	call FindActiveClusters(nActiveClusters, activeClusterIndices, particle, nClusters, iLastLocus) ! --in findfunctions.F90
	if (nActiveClusters > 0) then				
		do iActiveCluster = 1, nActiveClusters
			iCluster = activeClusterIndices(iActiveCluster)	
			if (particle(iCluster)%phase /= 1) cycle	
			call GetParticleAverageAttributes(nImagingDeployDepths, attsInSizeClass(:,:,:,iProfile), &
				attsInVeloClass(:,:,:,iProfile), attsInMainType(:,:,:,iProfile), particle, &
				nClusters, iCluster, iYear) ! --in modeloutput.F90						
#ifdef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS
			call GetParticleAverageAttributesByVolumeClass(attsInVolumeClass(:,:,iProfile), &
				particle, nClusters, iCluster, iTimeStep) ! --in modeloutput.F90	
#endif
		end do			
	end if

!#ifndef BLOCK_TMM_SLAMS
!#ifndef SAVE_DISK_SPACE
	! If this is the last model year, record particle attributes of all clusters. 
	if (iYear == nYears .and. isMidMonth .and. iLastLocus > 0) then
		call FindPhase1clusters(nPhase1clusters, phase1clusterIndices, particle, nClusters, iLastLocus) ! --in findfunctions.F90
		if (nPhase1clusters > 0) then
			call WriteSnapshotsOfParticleAttributes(nPhase1clusters, phase1clusterIndices(:), &
				particle, nClusters, nProfiles) ! --in modeloutput.F90	
		end if
		write(*,*)
 		write(*,*) 'Write mid-month snapshot'
 		write(*,*) 'Year:', iYear, ' Month:', m
 		write(*,*) 'iTimeStepYear:', iTimeStepYear
 		write(*,*)
		call WriteInstantaneousSnapshots(localSeafloorDepth,nProfiles) ! --in modeloutput.F90
	end if
!#endif                 
!#endif

#ifndef BLOCK_FEW_TIME_STEPS	
	end if ! doDiagnosticsAccum
#endif

#endif ! BLOCK_DIAGNOSTICS

! ----------------------------------------------------------------------------------------
! ------------------------------------DEPTH LOOP BEGINS-----------------------------------
! ----------------------------------------------------------------------------------------

	nZooDeadClustersPerProfile(iProfile) = 0
	call ZooplanktonDielVerticalMigrationDepth(dvmDepthNightUpper, dvmDepthNightLower, &
		dvmDepthDayUpper, dvmDepthDayLower, dvmBoundsBianchi(:), nSurfLayersBianchi, &
		nMesoLayersBianchi, SMSterm(:,:,iProfile), auxTerm(:,1,iProfile), auxCount(:,1,iProfile), &
		nLocalDepthLayers, ztop(:), zbot(:), O2(:,iForcingStep), TempC(:,iForcingStep), &
		Chla(iForcingStep), MLD(iForcingStep), localSeafloorDepth) ! --in heterotrophicmetabolism.F90

	do iDepthLayer = 1, nLocalDepthLayers
		zmidLayer = FindDepthLayerMidpointDepth(iDepthLayer, nLocalDepthLayers, ztop(:), zbot(:))
		ztopLayer = ztop(iDepthLayer) ! m
		zbotLayer = zbot(iDepthLayer) ! m
		gridCellVolume = (zbotLayer-ztopLayer)*gridCellArea(iProfile) ! m3
		tkeLayer = tkeDepthProfile(iDepthLayer)

		! Identify clusters in the current depth layer that have phase=1
		call FindClustersInDepthLayer(nLayerClusters, layerClusterIndices, particle, &
			nClusters, nActiveClusters, activeClusterIndices(:), ztopLayer, zbotLayer) ! --in findfunctions.F90

		! Process clusters in this depth layer
		if (nLayerClusters == 0) cycle ! skip if no clusters found
		call ParticleProcessing(particle, nClusters, iLastLocus, SMSterm(:,iDepthLayer,iProfile), &
			auxTerm(:,iDepthLayer,iProfile), auxCount(:,iDepthLayer,iProfile), nLayerClusters, &
			layerClusterIndices(:), PAR0(iForcingStep), OmegaCalc(iDepthLayer,iForcingStep), &
			O2(iDepthLayer,iForcingStep), TempC(iDepthLayer,iForcingStep), Rho(iDepthLayer,iForcingStep), &
			MesoZoo(iDepthLayer,iForcingStep), DynVisco(iDepthLayer,iForcingStep), tkeLayer, &
			dvmDepthNightUpper, dvmDepthNightLower, dvmDepthDayUpper, dvmDepthDayLower, &
			zmidLayer, gridCellVolume, gridCellArea(iProfile), iTimeStep, iYear, iProfile) ! --in wrappers.F90							
	end do

! ----------------------------------------------------------------------------------------
! ------------------------------------DEPTH LOOP ENDS-------------------------------------
! ----------------------------------------------------------------------------------------

#ifndef BLOCK_TMM_SLAMS
	! Update particle array as zooplankton dead body clusters may have been added
	! during the depth loop. Find clusters in particle array with id>0.
	call FindActiveClusters(nActiveClusters, activeClusterIndices, particle, nClusters, iLastLocus) ! --in findfunctions.F90
#endif

! ----------------------------------------------------------------------------------------
! -------------------------------SINKING-RELATED PROCESSES--------------------------------
! ----------------------------------------------------------------------------------------

	if (nActiveClusters > 0) then
		call GravitationalSettling(particle, nClusters, nActiveClusters, activeClusterIndices(:), &
			SMSterm(:,1,iProfile), nSedTrapDeployDepths, sedTrap(:,:,iProfile), &
			sedTrapSf(:,:,iProfile), sedTrapDeployDepths(:), zeu, MLD(iForcingStep), &
			localSeafloorDepth, nLocalDepthLayers, ztop(:), zbot(:), Rho(1,iForcingStep), &
			DynVisco(1,iForcingStep), muPftProfile(:,:), doDiagnosticsAccum, iTimeStep, iYear, iProfile) ! --in wrappers.F90

#ifdef BLOCK_DIAGNOSTICS
#ifndef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS			
		if (doDiagnosticsAccum) then
#endif
		! Process clusters that have reached the seafloor (phase=2)
		do iActiveCluster = 1, nActiveClusters
			iCluster = activeClusterIndices(iActiveCluster)		
			if (particle(iCluster)%phase /= 2) cycle		
			call GetParticleAverageAttributes(1, attsInSizeClassSf(:,:,:,iProfile), &
				attsInVeloClassSf(:,:,:,iProfile), attsInMainTypeSf(:,:,:,iProfile), &
				particle, nClusters, iCluster, iYear) ! --in modeloutput.F90	
			particle(iCluster)%phase = 3 ! after collection by the traps, release the cluster	
		end do
#ifndef BLOCK_FIXED_INITIAL_PARTICLE_CONDITIONS	
		end if ! doDiagnosticsAccum
#endif
#endif			
	end if ! nActiveClusters > 0

! ----------------------------------------------------------------------------------------
! ---------------UPDATE SMS, PARTICLE ARRAY AND WRITE PERIODIC DIAGNOSTICS----------------
! ----------------------------------------------------------------------------------------
		
	! SMS term update.
	call SourcesMinusSinksTermInformation(nLocalDepthLayers, iProfile) ! --in modeleulerianvariables.F90

#ifdef BLOCK_SQUEEZE_PARTICLE_ARRAY
#ifndef BLOCK_TMM_SLAMS

	! Squeeze the particle array to speed up the model. After nTimeStepsPassedToShiftClusters, 
	! shift clusters to occupy spaces that have been left by the clusters that reached 
	! the seafloor, or were dissolved (tiny clusters), or were emptied after bacteria 
	! respiration, or were emptied after mineral dissolution, and which have phase == 3.											
	if (MOD(iTimeStep,nTimeStepsPassedToShiftClusters) == 0 .or. iTimeStep == maxNumTimeSteps) then
		call FindPhase1clusters(nPhase1clusters, phase1clusterIndices, particle, nClusters, iLastLocus) ! --in findfunctions.F90
		call CompactParticleArray(particle, nClusters, iLastLocus, nPhase1clusters, phase1clusterIndices(:)) ! --in calcparticleattributes.F90
	end if
		
#endif
#endif
	else ! end checking whether there has been at least one cluster created at the surface in this time step
		write(*,*)
		write(*,'(a)') 'No clusters seeded: forcing conditions not met'
		write(*,'(a,es12.4,a)') '  NPP:        ', NPP(iForcingStep), ' mol C m-2 s-1'
		write(*,'(a,f6.2,a)')  '  Temperature:', TempC(1,iForcingStep), ' °C'
		write(*,'(a,es12.4,a)') '  PAR0:       ', PAR0(iForcingStep), ' W m-2'
		write(*,'(a,f7.2,a)')  '  NO3@surf:   ', NO3(1,iForcingStep), ' mmol m-3'
		write(*,*)
	end if ! iLastLocus > 0 / nExistingClusters + nNewClusters > 0
	
#ifdef BLOCK_DIAGNOSTICS
#ifndef BLOCK_TMM_SLAMS
#ifndef BLOCK_FEW_TIME_STEPS
	! Write periodic information to files/print to the screen	
	if (doDiagnosticsAccum) then
		if (iProfile == nProfiles) then
			diagTimer%count = diagTimer%count + 1 ! increment the count at every time step
			if (diagTimer%count==diagTimer%numTimeSteps) then ! after 1 month, time to write averages to file
				doDiagnosticsAvg = .TRUE.
				diagsAccumTimeSteps = diagTimer%count
			end if
#endif
#endif 

#ifdef BLOCK_FEW_TIME_STEPS	
			diagsAccumTimeSteps	= 1	
			call WritePeriodicOutput(particle, nClusters, iLastLocus, localSeafloorDepth, &
				diagsAccumTimeSteps, nProfiles, iTimeStep) ! --in wrappers.F90
			write(*,*)
			write(*,'(a,i0,a,i0,a)') 'Writing diagnostics at tstep ', iTimeStep, ', averaged over ', diagsAccumTimeSteps, ' tsteps'
			write(*,*)
#else
			if (doDiagnosticsAvg) then
				call WritePeriodicOutput(particle, nClusters, iLastLocus, localSeafloorDepth, &
					diagsAccumTimeSteps, nProfiles, iTimeStep) ! --in wrappers.F90
				write(*,*)
				write(*,'(a,i0,a,i0,a)') 'Writing diagnostics at tstep ', iTimeStep, ', averaged over ', diagsAccumTimeSteps, ' tsteps'
				write(*,*)
			end if
#endif           
					
#ifndef BLOCK_TMM_SLAMS
#ifndef BLOCK_FEW_TIME_STEPS
			! Reset the timer for the next averaging interval
			if (doDiagnosticsAvg) then
				call UpdateStepTimer(iTimeStep, diagTimer) 
				doDiagnosticsAvg = .FALSE.
			end if
		end if ! iProfile==nProfiles
	end if ! doDiagnosticsAccum
#endif
#endif	
#endif

end do ! end time step loop

! ----------------------------------------------------------------------------------------
! -------------------------------------TIME STEP ENDS-------------------------------------
! ----------------------------------------------------------------------------------------

call cpu_time( stopTimer )
write(*,*)
write(*,'(a,f0.2,a)') 'Model time is ', (stopTimer-startTimer)/60d0, ' minutes.'
	
#ifdef BLOCK_DIAGNOSTICS
#ifndef BLOCK_TMM_SLAMS
if (iLastLocus > 0) then
	call WriteModelClosureInformation(iTimeStep, iLastLocus) ! --in modeloutput.F90
	call WriteStatus('completed')
	write(*,*)
	write(*,*) '============================== THE END ==================================='
	write(*,*)
else
	call WriteStatus('no clusters seeded')
end if
#endif
#endif

end program BiologicalCarbonPump