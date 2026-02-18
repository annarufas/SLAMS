module timer

use modelparameters, only: nTimeStepsDay

implicit none

private
public :: StepTimer, IniStepTimer, UpdateStepTimer

integer, dimension(12) :: daysInMonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

type :: StepTimer
  logical :: fixedStep, haveResetStartTimeStep
  integer :: count, startTimeStep, startTimeStepResetFreq, numTimeSteps, maxNumIntervals, currInterval
  integer, dimension(12) :: timeIntervals
end type StepTimer

contains

! ========================================================================================

subroutine IniStepTimer(Iter0, theTimer)

	integer, intent(in) :: Iter0
	type(StepTimer), intent(inout) :: theTimer

    integer :: it 
!     integer, dimension(12) :: tmparr

!   Ultimately want to read this from a namelist
    thetimer%startTimeStep = Iter0 + 1; ! by default we start at first time step
!     ierr = PetscOptionsGetInt(NULL,pre,"-start_time_step",&thetimer%startTimeStep,&flg);CHKERRQ(ierr);
! 	ierr = PetscPrintf(PETSC_COMM_WORLD,"Start time step for StepTimer object %s is %d\n", pre, thetimer%startTimeStep);CHKERRQ(ierr);	  

    thetimer%maxNumIntervals = 12
!     ierr = PetscOptionsGetIntArray(NULL,pre,"-time_steps",tmparr,&thetimer%maxNumIntervals,&flg);CHKERRQ(ierr);
!     if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of step timer time steps with the -%stime_step flag",pre);

    if (thetimer%maxNumIntervals==1) then
      thetimer%fixedStep = .TRUE.
	  thetimer%currInterval = 1 ! Not used but we set it anyway to be safe
      thetimer%numTimeSteps = nTimeStepsDay*30 ! tmparr(1)
    else
      thetimer%fixedStep = .FALSE.      
! 	  PetscMalloc(thetimer%maxNumIntervals*sizeof(PetscInt), &thetimer%timeIntervals);
!       ierr = PetscPrintf(PETSC_COMM_WORLD,"Variable number of intervals specified for StepTimer object %s\n", pre);CHKERRQ(ierr);	  
	  thetimer%timeIntervals(1:thetimer%maxNumIntervals) = nTimeStepsDay*daysInMonth(1:thetimer%maxNumIntervals)
! 		ierr = PetscPrintf(PETSC_COMM_WORLD,"  Interval #%d=%d\n", it+1,thetimer%timeIntervals[it]);CHKERRQ(ierr);        
	  thetimer%currInterval = 1
	  thetimer%numTimeSteps = thetimer%timeIntervals(thetimer%currInterval)
    end if

    thetimer%startTimeStepResetFreq = -1
!     ierr = PetscOptionsGetInt(NULL,pre,"-start_time_step_reset_freq",&thetimer%startTimeStepResetFreq,&flg);CHKERRQ(ierr);
!     if (flg) {
! 	  PetscInt tmp=0;
! 	  if (!thetimer%fixedStep) {
! 		for (it=0; it<thetimer%maxNumIntervals; it++) {
! 		  tmp=tmp+(thetimer%timeIntervals[it]);
!         } 
! 	  } else {
!         tmp=thetimer%numTimeSteps;
!       }
!       if (tmp > thetimer%startTimeStepResetFreq) {
!         SETERRQ1(PETSC_COMM_WORLD,1,"Start time reset frequency less than total number of timer steps for StepTimer object %s",pre);
!       } 
!       ierr = PetscPrintf(PETSC_COMM_WORLD,"Start time will be reset every %d steps for StepTimer object %s\n", thetimer%startTimeStepResetFreq, pre);CHKERRQ(ierr);	  
!     }    
    
    thetimer%haveResetStartTimeStep = .FALSE.
	thetimer%count=0

end subroutine IniStepTimer

! ========================================================================================

subroutine UpdateStepTimer(Iter, theTimer)

	integer, intent(in) :: Iter
	type(StepTimer), intent(inout) :: theTimer

    logical :: endOfSequence = .TRUE.
    
!     ierr = PetscPrintf(PETSC_COMM_WORLD,"Updating StepTimer object %s at iter %d\n", pre,Iter);CHKERRQ(ierr);        

	thetimer%count = 0 ! reset counter
    thetimer%haveResetStartTimeStep = .FALSE.
    
    if (.NOT.(thetimer%fixedStep)) then
      thetimer%currInterval = thetimer%currInterval + 1
      if (thetimer%currInterval>thetimer%maxNumIntervals) then
!       We're now at the end of the sequence     
        thetimer%currInterval = 1
      else
!       Still within sequence      
        endOfSequence = .FALSE.
      end if       
      thetimer%numTimeSteps = thetimer%timeIntervals(thetimer%currInterval)
!       ierr = PetscPrintf(PETSC_COMM_WORLD,"New interval for StepTimer object %s at iter %d is %d\n", pre,Iter,thetimer%numTimeSteps);CHKERRQ(ierr);
    end if

    if ((thetimer%startTimeStepResetFreq > 0) .AND. (endOfSequence)) then
      thetimer%startTimeStep = thetimer%startTimeStep + thetimer%startTimeStepResetFreq
      thetimer%haveResetStartTimeStep = .TRUE.
!       ierr = PetscPrintf(PETSC_COMM_WORLD,"New start time step for StepTimer object %s is %d\n", pre, thetimer%startTimeStep);CHKERRQ(ierr);
    end if

end subroutine UpdateStepTimer
						
! ========================================================================================

end module timer