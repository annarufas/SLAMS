module modelcounters

! ----------------------------------------------------------------------------------------
! This module is concerned with initialising and carrying the various model counters that
! are used to calculate model diagnostics.
! ----------------------------------------------------------------------------------------

implicit none

public
					
integer                              :: iYear,              &
			 							iTimeStepYear

real*8                               :: particleCollectionPeriod ! days

integer                              ::	iLastSliceAvgAtt,    &
			 							iLastSliceAvgVolAtt, &
			 							iLastSliceSedTrap,   &
			 							iLastSliceEulerian,  &
			 							iLastSliceAnnual,    &
			 							iLastSliceAux,       &
			 							iLastSliceStats,     &
			 							iLastSliceSnapshot,  &
			 							iLastSliceParticleNumSnapshot

integer*8, dimension(:), allocatable :: nTimesEnteringCollLoop,          &
			 							nCollisionsType1,                &
			 							nCollisionsType2,                &
			 							nCollisionsType3,                &
			 							nCollisionsType4,                &
			 							nParticlesEvaluatedForGrazing,   &
			 							nZooEncounteredParticles,	     &
			 							nZooIngestionEvents,             &
			 							nZooFragmentationEvents,         &
			 							nZooNegletedFoodParticles,       &
			 							nFaecalPelletClustersProduced,   &
										nZooDeadClustersProduced,        &
			 							nZooDeadClustersPerProfile,      &
			 							nMicrobRespiredClusters,         &              
			 							nFragmentedBigClusters,          &
			 							nRespiredTinyClusters,           &
			 							nMineralClustersDissolved,       &
			 							nClustersPhotolysed,             &
			 							nClustersAtSeafloor	 							                      
		
contains

! ========================================================================================	

subroutine InitialiseCounters( nProfiles )

	integer, intent(in) :: nProfiles
		 								
	allocate(nTimesEnteringCollLoop(nProfiles))
	allocate(nCollisionsType1(nProfiles))
	allocate(nCollisionsType2(nProfiles))
	allocate(nCollisionsType3(nProfiles))
	allocate(nCollisionsType4(nProfiles))
	allocate(nParticlesEvaluatedForGrazing(nProfiles))
	allocate(nZooEncounteredParticles(nProfiles))
	allocate(nZooIngestionEvents(nProfiles))
	allocate(nZooFragmentationEvents(nProfiles))
	allocate(nZooNegletedFoodParticles(nProfiles))
	allocate(nFaecalPelletClustersProduced(nProfiles))
	allocate(nZooDeadClustersProduced(nProfiles))
	allocate(nZooDeadClustersPerProfile(nProfiles))
	allocate(nMicrobRespiredClusters(nProfiles))
	allocate(nFragmentedBigClusters(nProfiles))
	allocate(nRespiredTinyClusters(nProfiles))
	allocate(nMineralClustersDissolved(nProfiles))
	allocate(nClustersPhotolysed(nProfiles))
	allocate(nClustersAtSeafloor(nProfiles))
	
	particleCollectionPeriod = 0d0
	iLastSliceAvgAtt    = 0
	iLastSliceAvgVolAtt = 0
	iLastSliceSedTrap   = 0
	iLastSliceEulerian  = 0
	iLastSliceAnnual    = 0
	iLastSliceAux       = 0
	iLastSliceStats     = 0
	iLastSliceSnapshot  = 0
	iLastSliceParticleNumSnapshot = 0
			
	iYear              = 0
	iTimeStepYear      = 0

	nTimesEnteringCollLoop(:)          = 0
	nCollisionsType1(:)                = 0
	nCollisionsType2(:)                = 0
	nCollisionsType3(:)                = 0
	nCollisionsType4(:)                = 0
	nParticlesEvaluatedForGrazing(:)   = 0
	nZooEncounteredParticles(:)        = 0
	nZooIngestionEvents(:)             = 0
	nZooFragmentationEvents(:)         = 0
	nZooNegletedFoodParticles(:)	   = 0
	nFaecalPelletClustersProduced(:)   = 0
	nZooDeadClustersProduced(:)        = 0
	nZooDeadClustersPerProfile(:)      = 0
	nMicrobRespiredClusters(:)         = 0
	nFragmentedBigClusters(:)          = 0
	nRespiredTinyClusters(:)           = 0
	nMineralClustersDissolved(:)       = 0
	nClustersPhotolysed(:)             = 0
	nClustersAtSeafloor(:)             = 0
			 
end subroutine InitialiseCounters

! ========================================================================================	

end module modelcounters