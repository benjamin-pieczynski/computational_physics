! Program: schrodinger
! By: Benjamin Pieczynski
!-----------------------------------------------------------------------------
! This program determines the solutions for the time-dependent schrodinger
! equation for either the case of harmonice or infinite well (v=0). The 
! program uses the Crank-Nicolson method for the determination of the 
! Schrodinger equation. This method will split the discretized hamiltonian
! into half implict and half explicit and uses linear algebra to solve. More
! information can be found within the readme.
!-----------------------------------------------------------------------------
program schrodinger 

use types
use read_write, only : read_input, write_time_evolution, write_expectation_values
use quantum, only : sample_box, construct_initial_wavefunction, construct_time_evolution_matrix, &
    evolve_wave_function, expectation_values

implicit none

real(dp) :: length, delta_t, width, center, k_oscillator
integer :: n_points, n_steps
character(len=1024) :: time_file, density_file 
real(dp), allocatable :: x_vector(:) !will be of size n_points.
real(dp), allocatable :: wave_func(:)! will be of size 2*n_points.
real(dp), allocatable :: evolution_matrix(:,:) !will be of size 2*n_points by 2*n_points
real(dp), allocatable :: time_wave_func(:,:) !will be of size n_points by n_steps + 1 (the +1 is so that you can store the t=0 value)
real(dp), allocatable :: norm(:), position(:), sigma(:) !all of size n_steps + 1 

call read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator, time_file, density_file)
call sample_box(length, n_points, x_vector)
call construct_initial_wavefunction(n_points, center, width, wave_func, x_vector)
call construct_time_evolution_matrix(n_points, length, k_oscillator, x_vector, delta_t, evolution_matrix)
call evolve_wave_function(n_points, n_steps, wave_func, evolution_matrix, time_wave_func)
call expectation_values(n_points, length, n_steps, x_vector, time_wave_func, norm, position, sigma)
call write_time_evolution(n_steps, x_vector, time_wave_func, density_file)
call write_expectation_values(n_steps, delta_t, norm, position, sigma, time_file)

end program schrodinger