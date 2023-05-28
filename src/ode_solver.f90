!-----------------------------------------------------------------------
!Module: ode_solver
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module uses a 4th order Runge Kutta method to solve for the 
!! 8 coupled ordinary differential equations specified in planets_ode.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_runge_kutta_4
!!----------------------------------------------------------------------
module ode_solver
use types
implicit none
private

public :: solve_runge_kutta_4

interface
    function func(r, t, work) result(f)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: r(:), t, work(:)
        real(dp), allocatable :: f(:)
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Subroutine: solve_runge_kutta_4
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine solves the 8 coupled ordinary differential equations
!! specified in the system. They are solved using a 4th order Runge Kutta
!! method. The output is the resulting position and velocities over time.
!!----------------------------------------------------------------------
!! Input:
!!
!! f        procedure       function containing ODEs to solve
!! t_f      real            final time
!! r_init   real            initial conditions
!! n        integer         step size
!! work     real            array containing object masses
!!----------------------------------------------------------------------
!! Output:
!!
!! r        real            solution array
!! t        real            time array
!-----------------------------------------------------------------------
subroutine solve_runge_kutta_4(f, t_f, r_init, n, work, r, t)
    implicit none
    procedure(func) :: f
    integer, intent(in) :: n
    real(dp), intent(in) :: t_f, r_init(:)
    real(dp), allocatable, intent(out) :: t(:), r(:,:)
    real(dp), intent(in) :: work(:)
    real(dp), allocatable :: r_solve(:), k1(:), k2(:), k3(:), k4(:)
    real(dp) :: h, t_i, t_s
    integer :: i, num_var

    ! allocating memory
    num_var = size(r_init)
    allocate(r_solve(1:n))
    allocate(t(1:n))
    allocate(r(1:num_var, 1:n))
    allocate(k1(1:n))
    allocate(k2(1:n))
    allocate(k3(1:n))
    allocate(k4(1:n))

    ! initializing variables
    t_i = 0._dp
    h = (t_f - t_i)/n
    r_solve = r_init
    t_s = t_i
    
    do i=1,n
        k1 = h*f(r_solve, t_s, work)
        k2 = h*f(r_solve + k1/2, t_s + h/2, work)
        k3 = h*f(r_solve + k2/2, t_s + h/2, work)
        k4 = h*f(r_solve + k3, t_s + h, work)
        r_solve = r_solve + (k1 + 2*k2 + 2*k3 + k4)/6
        t_s = t_s + h
        r(:,i) = r_solve
        t(i) = t_s
    end do

end subroutine solve_runge_kutta_4

end module ode_solver