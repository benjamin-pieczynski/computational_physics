!-----------------------------------------------------------------------
!Module: analytic_functions
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module creates the functions that will be called for the
!! analytic solution for x*sin(x) and its second derivative. The 
!! functions will use the real user input x_zero in order to calculate
!! the analytic solution.  The analytic solution will be the output
!! y_zero.  The euler_formulas.f90 module will utilize the analytic_f(x)
!! function module to approximate the second derivative of x*sin(x).
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!! analytic_f(x)
!! second_derivative_f(x)
!-----------------------------------------------------------------------
module analytic_functions
use types
implicit none

! The private statement restricts every function, parameter and variable
! defined in this module to be visible only by this module
private
! Then we use the public statement to only make visible to other modules 
! the few functions or subroutines that will be used by them
public analytic_f, second_derivative_f

contains

!-----------------------------------------------------------------------
!Function: analytic_f
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function takes in and evaluates the point x_zero in the
!! mathematical expression f(x_0) = x_0 sin(x_0).  The result y_zero
!! is the analytic solution for f(x_0).  This value will be used in
!! euler_formulas.f90 to solve for the 3 and 5 point approximations
!! of f''(x_0).
!!
!!----------------------------------------------------------------------
!! Arguments:
!!
!! x_zero	real	point x_0 at which to evaluate f(x_0)
!-----------------------------------------------------------------------
!! Result:
!!
!! y_zero	real	x_0 sin(x_0)
!-----------------------------------------------------------------------
function analytic_f(x_zero) result(y_zero)
    implicit none
    real(dp), intent(in) :: x_zero
    real(dp) :: y_zero
    ! The function should return
    ! x*sin(x)
    y_zero = x_zero*sin(x_zero)
end function analytic_f

!-----------------------------------------------------------------------
!Function: second_derivative_f
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function takes in and evaluates the point x_zero for the
!! mathematical expression f''(x_0) = 2cos(x_0) - x_0 sin(x_0).  The result
!! y_zero serves as the analytic solution for f''(x_0).  This value will
!! be compared to the euler approximations for f''(x_0).
!!
!!----------------------------------------------------------------------
!! Arguments:
!!
!! x_zero	real	point x_0 at which to evaluate f''(x_0)
!-----------------------------------------------------------------------
!! Result:
!!
!! y_zero	real	
!-----------------------------------------------------------------------
function second_derivative_f(x_zero) result(y_zero)
    implicit none
    real(dp), intent(in) :: x_zero
    real(dp) :: y_zero
    ! the second derivative of x*sin(x) and stored in y_zero
    y_zero = (2*cos(x_zero) - x_zero*sin(x_zero))
end function second_derivative_f
    
end module analytic_functions
