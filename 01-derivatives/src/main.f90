! Program: numerical_derivatives
! By: Benjamin Pieczynski
!-----------------------------------------------------------------------------
!
! Calculates the second derivative of f(x) = x sin(x) using Euler's method.
! Euler's method uses the 3 point and 5 point methods for approximating the
! second derivative:
!
! 5 point method: f''(x_0) ~ (-f(x-2h)+16f(x-h)-30f(x)+16f(x+h)-f(x+2h))/12h^2
! 3 point method: f''(x_0) ~ (f(x+h)-2f(x)+f(x-h))/(h^2)
!
! The program will prompt the user for a real number input for x_zero. Then
! the program will solve for the analytic solution and approximate the 3 and
! 5 point Euler methods. The results will be exported to results.dat.  To
! accomplish this the program uses four seperate modules. 
! 
! (1) The program calls read_input from the read_write.f90 module. This will 
! ask for a user input(x_zero) and check to see if the input is a real number. 
! 
! (2) The program will call the write_derivatives subroutine in that same
! module. This subroutine will calculate h_step and utilize a do loop to
! iterate until the h_step * h_increment values are greater than h_max.
! The subroutine also solves the analytic solution found in module
! analytic_functions.f90.
!
! (3) Inside the loop, the module calls the euler_5points and euler_3points
! functions found in euler_formulas.f90. These functions will calculate the
! approximations for the second derivatives. After the values are calculated
! they are written to to the results.dat file befor the loop reiterates.
!
! (4) Both the euler_3points and euler_5points use module 
! analytic_functions.f90's function analytic_f to find f(x) for the 3 an 5
! point Euler methods.
!
! (5) The subroutine exits the loop when h_step > h_max. The subroutine and
! module end, thus causing the program to end with the output of data sent to
! results.dat.
!
!-----------------------------------------------------------------------------
program numerical_derivatives
use types 
use read_write

implicit none
real(dp) :: x_zero

! the dp inside of the real(dp) declaration was defined in the types module
! and allows to use reals with double precision 
call read_input(x_zero)
call write_derivatives(x_zero)
end program numerical_derivatives