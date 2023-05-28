!-----------------------------------------------------------------------------
!! Program: schrodinger_solution
!! By: Benjamin Pieczynski
!!
!! This program determines the numerical and analytical solutions to the one
!! dimesional schrodinger equation over different potentials. The program uses
!! a lapack library to solve for the energies(eigenvalues) and the wavefunctions
!! (eigenvectors). Solutions are produced for the Infinite Well, Harmonic
!! Oscillator, and the Woods-Saxon. The user inputs a length, number of iterations,
!! and radius for the woods saxon potential. The results are written to different
!! files for each potential. The program ultimately finds the energies and 
!! wavefunctions for the Woods-Saxon potential.
!!
!! particle in a box: infinite_well_wf.dat
!! Harmonic oscillator: harmonic_oscillator_wf.dat
!! Woods-Saxon: woods_saxon_wf.dat
!! Woods-Saxon Energies: woods_saxon_ener.dat
!!
!-----------------------------------------------------------------------------
program schrodinger_solution 

use types
use read_write, only : read_input, write_probability_density, print_energies, write_woods_saxon_energies
use qm_solver, only: sample_box, solve_infinite_well, analytic_infinite_well, solve_harmonic_oscillator,&
    analytic_harmonic_oscillator, solve_woods_saxon
implicit none

integer :: n_points
real(dp) :: length
integer :: radius

integer, parameter :: n_energies = 3
real(dp) :: energies(1:n_energies), analytic_energies(1:n_energies)
real(dp), allocatable :: wave_functions(:,:)
real(dp), allocatable :: x_vector(:)
integer, parameter :: r_min = 2, r_max = 10

call read_input(n_points, length, radius)
call sample_box(length, n_points, x_vector)

! Solving particle in a box
call solve_infinite_well(n_points, length, energies, wave_functions, x_vector)
call analytic_infinite_well(length, analytic_energies)
call print_energies('Infinite Well', energies, analytic_energies)
call write_probability_density('infinite_well_wf.dat', x_vector, wave_functions)

! Solving harmonic oscillator
call solve_harmonic_oscillator(n_points, length, energies, wave_functions, x_vector)
call analytic_harmonic_oscillator(analytic_energies)
call print_energies('Harmonic oscillator', energies, analytic_energies)
call write_probability_density('harmonic_oscillator_wf.dat', x_vector, wave_functions)

! Solving Woods Saxon
call solve_woods_saxon(n_points, length, radius, energies, wave_functions, x_vector)
call write_probability_density('woods_saxon_wf.dat', x_vector, wave_functions)

! Woods Saxon Energies as a function of radius
call write_woods_saxon_energies('woods_saxon_ener.dat', n_points, length, r_min, r_max, x_vector)

end program schrodinger_solution