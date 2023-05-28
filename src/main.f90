! Program: planets
! By: Benjamin Pieczynski
!-----------------------------------------------------------------------------
! This program can model the orbits of two planets or the orbits of a planet
! and moon about a primary object (star). The program solves 8 coupled first
! order ordinary differential equtions using a 4th-order runge-kutta method.
! The program has a set of default parameters, but passing a namelist into the
! program allows a user to test different initial conditions for the system.
!-----------------------------------------------------------------------------
program planets

use types
use read_write, only : read_input, write_results
use ode_solver, only : solve_runge_kutta_4
use mechanics, only : calculate_energy, calculate_iomega, planets_ode

implicit none

real(dp) :: work_array(1:3), initial_condition(1:8)
real(dp) :: final_time
integer :: n_steps
real(dp), allocatable :: time(:), solution(:,:), energy(:), ang_momentum(:)
character(len=1024) :: output_file

call read_input(work_array, initial_condition, final_time, n_steps, output_file)

call solve_runge_kutta_4(planets_ode, final_time, initial_condition, n_steps, work_array, solution, time)
call calculate_energy(solution, n_steps, work_array, energy)
call calculate_iomega(solution, n_steps, work_array, ang_momentum)

call write_results(time, solution, energy, ang_momentum, n_steps, output_file)


end program planets