module montecarlosampling

! ----------------------------------------------------------------------------------------
! This module is concerned with the Monte Carlo Method.
! ----------------------------------------------------------------------------------------

use, intrinsic :: ieee_arithmetic
use sanitychecks, only: WriteStatusAndStop	

implicit none
private
public :: RandProbCase

contains

! ========================================================================================

function RandProbCase( nCases, listOfCases )

	!-------------------------------------------------------------------------------------
	! This function returns the location of an object randomly picked up from a cumulative 
	! distribution function (CDF) of the different possible objects.
	!-------------------------------------------------------------------------------------

	integer :: RandProbCase
	integer, intent(in) :: nCases
	real*8, dimension(nCases), intent(in) :: listOfCases
	
	integer :: nCasesWithValues, iCase, nEdges, idx
	real*8 :: maxCumSumValue, harvest
	real*8, dimension(:), allocatable :: valuesInNonZeroCases, cumSumValues, edges	
	logical, dimension(:), allocatable :: matchCriteria 
	integer, dimension(:), allocatable :: casesWithValues

	if (SUM(listOfCases(:)) == 0d0) then
    	write(*,*) 'ERROR: in RandProbCase, all cases have zero probability'
    	call WriteStatusAndStop('invalid probabilities')
	end if

	! Initialise
	RandProbCase = 0 ! number that corresponds to the randomly drawn position on listOfCases
	allocate(matchCriteria(nCases))
	matchCriteria(:) = listOfCases(1:nCases) > 0
	nCasesWithValues = COUNT(matchCriteria(:))
	
	! Pack non-zero values in the vector (these are positions in the vector that could be empty)			
	if (nCasesWithValues < nCases) then ! if some positions are empty
		allocate(valuesInNonZeroCases(nCasesWithValues))
		valuesInNonZeroCases = PACK(listOfCases(:), matchCriteria(:))
	else if (nCasesWithValues == nCases) then
		allocate(valuesInNonZeroCases(nCases)) ! nCases = nCasesWithValues
		valuesInNonZeroCases(:) = listOfCases(:)
	end if

	! Compute cumulative sum and scale to [0, 1]
	allocate(cumSumValues(nCasesWithValues))
	cumSumValues(:) = 0d0
	cumSumValues(1) = valuesInNonZeroCases(1)
	do iCase = 2, nCasesWithValues
		cumSumValues(iCase) = cumSumValues(iCase-1) + valuesInNonZeroCases(iCase) 
	end do
	maxCumSumValue = cumSumValues(nCasesWithValues)
	cumSumValues(:) = cumSumValues(:) / maxCumSumValue ! SCALE

	! Construct edges for CDF
	allocate(edges(nCasesWithValues+1))	
	nEdges = nCasesWithValues + 1
	edges(1) = 0d0
	edges(2:nEdges) = cumSumValues(:)

	! Randomly pick up your object according to the cumulative distribution function	
	call random_number( harvest ) ! generates numbers in [0,1), not including 1
	do iCase = 1, nCasesWithValues			
		if (harvest >= edges(iCase) .and. harvest <= edges(iCase+1)) then
			RandProbCase = iCase
			exit
		end if
	end do	
	
	! Map back to original, non-packed vector	
	if (nCasesWithValues /= nCases) then	
		allocate(casesWithValues(nCases))
		idx = 1
		do iCase = 1, nCases
			if (listOfCases(iCase) > 0) then
				casesWithValues(idx) = iCase
				idx = idx + 1
			end if
		end do		
		RandProbCase = casesWithValues(RandProbCase)
		deallocate(casesWithValues)		
	end if

	! Cleanup
	deallocate(matchCriteria,valuesInNonZeroCases,cumSumValues,edges)	

end function RandProbCase       

! ========================================================================================

function RandCase(nCases)

	integer :: RandCase
	integer, intent(in) :: nCases
	
	real*8, dimension(:), allocatable :: values, cumSumValues, cumSumScaledValues
	real*8, dimension(:), allocatable :: edges
	integer :: iCase, nEdges
	real*8 :: harvest, maxCumSumValue

	! Initialise
	RandCase = 0
	allocate(values(nCases))
	allocate(cumSumValues(nCases))		
	allocate(cumSumScaledValues(nCases))
	allocate(edges(nCases+1))	
	
	! Construct a vector with the scaled cumsum values	
	values = 1d0	
	cumSumValues(1) = values(1)
	do iCase = 2, nCases
		cumSumValues(iCase) = cumSumValues(iCase-1) + values(iCase) 
	end do
	
	! Reescale to have the values between 0 and 1 (probability distribution function)
	maxCumSumValue = cumSumValues(nCases)
	cumSumScaledValues(1:nCases) = cumSumValues(:) / maxCumSumValue

	nEdges = nCases + 1
	edges(1) = 0d0
	edges(2:nEdges) = cumSumScaledValues(:)

	! Randomly pick up your object according to the probability distribution function	
	call random_number( harvest ) ! generates numbers in [0,1), not including 1
	do iCase = 1, nCases
		if (harvest >= edges(iCase) .and. harvest < edges(iCase+1)) then
			RandCase = iCase
			exit
		end if
	end do

	! Cleanup	
	deallocate(values,cumSumValues,cumSumScaledValues,edges)

end function RandCase

! ========================================================================================

end module montecarlosampling