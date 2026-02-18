module safemath

implicit none
private
public :: SafeFloorNonNegative, SafeCeilingNonNegative, SafeMOD, SafeDivide, SafeLog, &
	SafeSqrt, SafeNINTNonNegative, SafeSplitIntegerLike

contains

! ========================================================================================

function SafeFloorNonNegative(x) result(y)

	real(8), intent(in) :: x
    real(8) :: y

    ! Detect NaN safely (works even for signaling NaN)
    if (x /= x) then
        error stop 'ERROR: NaN in SafeFloorNonNegative'
    end if

    ! Reject infinities explicitly
    if (abs(x) > huge(1d0)) then
        error stop 'ERROR: Inf or overflow in SafeFloorNonNegative:'
    end if

    if (x <= 0d0) then
        y = 0d0
    else
        y = DINT(x) ! AINT truncates toward zero; for x>0 this == FLOOR(x)
    end if
    
end function SafeFloorNonNegative

! ========================================================================================

function SafeCeilingNonNegative(x) result(y)

    real(8), intent(in) :: x
    real(8) :: y

    if (x /= x) then
        error stop 'ERROR: NaN in SafeCeilingNonNegative'
    end if

    if (abs(x) > huge(1d0)) then
        error stop 'ERROR: Inf or overflow in SafeCeilingNonNegative:'
    end if

    if (x <= 0d0) then
        y = 0d0
    else
        y = DINT(x)
        if (y < x) y = y + 1d0
    end if
    
end function SafeCeilingNonNegative

! ========================================================================================

function SafeMOD(x, y) result(z)

    real(8), intent(in) :: x, y
    real(8) :: z

    if (x /= x .or. y /= y) then
        error stop 'ERROR: NaN in SafeMOD'
    end if

    if (y == 0d0) then
        error stop 'ERROR: MOD by zero in SafeMOD:'
    end if

    if (abs(x) > huge(1d0) .or. abs(y) > huge(1d0)) then
        error stop 'ERROR: Inf in SafeMOD:'
    end if

    ! MOD(x,y) = x - y*floor(x/y); use AINT for safety
    z = DMOD(x, y)
    
end function SafeMOD

! ========================================================================================

function SafeDivide(numer, denom, default) result(z)

    real(8), intent(in) :: numer, denom, default
    real(8) :: z
    real(8) :: abs_num, abs_den

    ! Detect NaN safely (works for quiet & signaling NaN)
    if (numer /= numer .or. denom /= denom) then
        error stop 'ERROR: NaN in SafeDivide'
    end if

    ! Division by zero
    if (denom == 0d0) then
        error stop 'ERROR: division by zero in SafeDivide'
    end if

    abs_num = abs(numer)
    abs_den = abs(denom)

    ! Detect overflow *before* division
    if (abs_den < 1d0 .and. abs_num > abs_den * huge(1d0)) then
        error stop 'ERROR: overflow in SafeDivide'
    end if

    z = numer / denom
    
end function SafeDivide

! ========================================================================================

function SafeLog(x) result(y)

    real(8), intent(in) :: x
    real(8) :: y

    if (x /= x) then
        error stop 'ERROR: NaN in SafeLog'
    end if

    if (x <= 0d0) then
        error stop 'ERROR: LOG of non-positive number'
    end if

    if (abs(x) > huge(1d0)) then
        error stop 'ERROR: Inf in SafeLog'
    end if

    y = LOG(x)
    
end function SafeLog

! ========================================================================================

function SafeSqrt(x) result(y)

    real(8), intent(in) :: x
    real(8) :: y

    if (x /= x) then
        error stop 'ERROR: NaN in SafeSqrt'
    end if

    if (x < 0d0) then
        error stop 'ERROR: SQRT of negative number'
    end if

    if (abs(x) > huge(1d0)) then
        error stop 'ERROR: Inf in SafeSqrt'
    end if

    y = SQRT(x)
    
end function SafeSqrt

! ========================================================================================

function SafeNINTNonNegative(x) result(y)

    real(8), intent(in) :: x
    real(8) :: y

    if (x /= x) then
        error stop 'ERROR: NaN in SafeNINTNonNegative'
    end if

    if (abs(x) > huge(1d0)) then
        error stop 'ERROR: Inf in SafeNINTNonNegative'
    end if

    if (x <= 0d0) then
        y = 0d0
    else
        ! Nearest integer in REAL arithmetic
        y = DINT(x + 0.5d0)
    end if
    
end function SafeNINTNonNegative

! ========================================================================================

subroutine SafeSplitIntegerLike(N, A, B)

    real(8), intent(in)  :: N
    real(8), intent(out) :: A, B
    real(8) :: half

    ! Detect NaN or infinite
    if (N /= N .or. abs(N) > huge(1d0)) then
        error stop 'ERROR: NaN or Inf in SafeSplitIntegerLike'
    end if

    ! Detect integer-like
    if (abs(N - DINT(N)) > 1d-10) then
        write(*,*) 'ERROR: non-integer-like value in SafeSplitIntegerLike'
        write(*,*) '  N =', N
        error stop
    end if

    ! Sanity: size
    if (N < 2d0) then
        write(*,*) 'ERROR: cannot split N < 2 in SafeSplitIntegerLike'
        write(*,*) '  N =', N
        error stop
    end if

    ! ---------- Split ----------
    half = DINT(N / 2d0)   ! floor(N/2)
    B    = half
    A    = N - B           ! ensures A + B = N

    ! ---------- Postconditions ----------
    if (A < 1d0 .or. B < 1d0) then
        write(*,*) 'ERROR: invalid split in SafeSplitIntegerLike'
        write(*,*) '  A, B =', A, B
        error stop
    end if

    if (abs((A + B) - N) > 1d-10) then
        write(*,*) 'ERROR: split does not conserve N'
        write(*,*) '  N, A+B =', N, A+B
        error stop
    end if

end subroutine SafeSplitIntegerLike

! ========================================================================================

end module safemath