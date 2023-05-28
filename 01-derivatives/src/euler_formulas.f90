!-----------------------------------------------------------------------
!Module: euler_formulas
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module approximates the second derivative of x*sin(x) using both
!! the 5 point and 3 point Euler methods.  These functions use
!! numerical approximations to approximate the second derivative. The
!! module uses the h_step as an input, as determined from the
!! read_write.f90 module. This module is called called from the
!! write_derivatives subroutine in read_write.f90 . This module uses the
!! analytic_f function in analytic_functions.f90 in order to generate an
!! approximation for f''(x_0).  This occurs in both the euler3points and
!! euler_5points functions.
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!! euler_3points(x,h)
!! euler_5points(x,h)
!-----------------------------------------------------------------------
module euler_formulas
use types
use analytic_functions, only : analytic_f
implicit none

! The private statement restricts every function, parameter and variable
! defined in this module to be visible only by this module
private
! Then we use the public statement to only make visible to other modules 
! the few functions or subroutines that will be used by them
public euler_3points
public euler_5points

contains

!-----------------------------------------------------------------------
!Function: euler_3points
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function utilizes the Euler 3 point numerical approximation in
!! order to approximate the second derivative of x sin(x).  It uses the
!! inputs of x_zero and h_step to approximate the result y_zero. y_zero
!! is approximately equal to f''(x_0) ~ (f(x+h)-2f(x)+f(x-h))/(h^2) .
!! The function analytic_f in the module analytic_functions.f90 is used
!! to solve for y_zero.
!!
!!----------------------------------------------------------------------
!! Arguments:
!!
!! x_zero   real    point x_0 at which to evaluate f''(x_0)
!! h_step   real    step size in the numerical expression
!-----------------------------------------------------------------------
!! Result:
!!
!! y_zero   real    (f(x+h)-2f(x)+f(x-h))/(h^2)
!-----------------------------------------------------------------------
function euler_3points(x_zero,h_step) result(y_zero)
    implicit none
    real(dp), intent(in) :: x_zero, h_step
    real(dp) :: y_zero
    real(dp) :: f_plus, f_zero, f_minus
    ! This evaluates the analytic function defined in the analytic_functions
    f_plus = analytic_f(x_zero + h_step)
    f_zero = analytic_f(x_zero)
    f_minus = analytic_f(x_zero - h_step)

    ! Here you can use the evaluated values to calculate the numerical
    ! approximation to the second derivative
    y_zero = (f_plus - (2 * f_zero) + f_minus) / (h_step**2)
end function euler_3points

!-----------------------------------------------------------------------
!Function: euler_5points
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function utilizes the Euler 5 point numerical approximation in
!! order to approximate the second derivative of x sin(x).  It uses the
!! inputs of x_zero and h_step to approximate the result y_zero. y_zero
!! is approximately equal to f''(x_0) ~ 
!! (-f(x-2h)+16f(x-h)-30f(x)+16f(x+h)-f(x+2h))/12h^2 . The function 
!! analytic_f in the module analytic_functions.f90 is used to solve for 
!! y_zero.
!!
!!----------------------------------------------------------------------
!! Arguments:
!!
!! x_zero   real    point x_0 at which to evaluate f''(x_0)
!! h_step   real    step size in the numerical expression
!-----------------------------------------------------------------------
!! Result:
!!
!! y_zero   real    (-f(x-2h)+16f(x-h)-30f(x)+16f(x+h)-f(x+2h))/12h^2
!-----------------------------------------------------------------------
function euler_5points(x_zero, h_step) result(y_zero)
    implicit none
    real(dp), intent(in) :: x_zero, h_step
    real(dp) :: y_zero
    real(dp) :: f_stage1, f_stage2, f_plus, f_zero, f_minus
    ! This evaluates the analytic function defined in the analytic_functions
    f_stage1 = analytic_f(x_zero - (2 * h_step))
    f_stage2 = analytic_f(x_zero + (2 * h_step))
    f_plus = analytic_f(x_zero + h_step)
    f_zero = analytic_f(x_zero)
    f_minus = analytic_f(x_zero - h_step)
    
    ! Here you can use the evaluated values to calculate the numerical
    ! approximation to the second derivative
    y_zero = (-f_stage1 + (16 * f_minus) - (30 * f_zero) + (16 * f_plus) - f_stage2) / (12 * h_step**2)
end function euler_5points

end module euler_formulas