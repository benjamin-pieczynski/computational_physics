! Program: nuclear_energies
! By: Benjamin Pieczynski
!-----------------------------------------------------------------------------
! This program reads in all up to date values of measured binding energies
! from a user input EXPERIMENT_AME2016.dat in order to estimate the position of 
! the nuclear dripline and determine the neutron and proton separation energies. 
! Estimation is done by using the semi-empirical mass formula in order to estimate 
! the drip-lines, which bounds the regions with positive separation energies.
! To determine the best fit for the linear model, matrices are used 
! to solve the set of linear equations. The uncertainty is determined using the
! a covariance matrix to propogate the experimental uncertainty to any quantity
! calculated with the model parameters. This program with then use the theoretical 
! model to calculate the binding energy for any isotope determined by different values 
! for Z and N. This will include finding the positions of the positions of valley of 
! stability for the neutron dripline. The experimental values and uncertainties
! are written to a file named results.dat along with the theoretical values
! and uncertainties.
!-----------------------------------------------------------------------------
program nuclear_energies
use types
use read_write, only : read_exp_data, write_predictions, write_advanced
use nuclear_model, only : find_best_parameters, valley_of_stability, theoretical_driplines
implicit none

integer, allocatable :: n_protons(:), n_neutrons(:), N_STABLE(:), N_DRIP(:), protons(:)
real(dp), allocatable :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)

! BASIC
call read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
call find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
call write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)

! ADVANCED
call valley_of_stability(c_parameters, protons, N_STABLE)
call theoretical_driplines(c_parameters, N_DRIP)
call write_advanced(c_parameters)

!------------------------------------------------------------

end program nuclear_energies