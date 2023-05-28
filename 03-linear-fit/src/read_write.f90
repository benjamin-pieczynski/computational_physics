!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! The purpose of this module is to read in the experimental data from
!! the EXPERIMENT_AME2016.d file (in the function read_exp_data), and to 
!! write the variables and call other modules that are needed for the 
!! write_predictions subroutine and the advanced subroutine.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_exp_data
!! write_predictions
!! write_advanced
!!----------------------------------------------------------------------
module read_write

use types
use nuclear_model, only : semi_empirical_mass, semi_empirical_error, valley_of_stability, theoretical_driplines

implicit none

private
public :: read_exp_data, write_predictions, write_advanced

contains

!-----------------------------------------------------------------------
!! Subroutine: read_exp_data
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine reads in the experimental data from the file specified
!! by the user. In the case of an invalid file entry to the user that the 
!! file specified could not be found. The data that will be read by this 
!! subroutine is the number of protons  
!!----------------------------------------------------------------------
!! Output:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
subroutine read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
    implicit none
    integer, intent(out), allocatable :: n_protons(:), n_neutrons(:)
    real(dp), intent(out), allocatable :: exp_values(:), uncertainties(:)
    character(len=128) :: filename, Nuc, A, N_exp, EA
    logical :: file_exists
    integer :: file_unit, n_points, i, skip_lines, stat, i_do
    real(dp) :: exp_values_r, uncertainties_r
    integer :: n_protons_r, n_neutrons_r


    ! use the print statement to display a message explaining what the program
    ! does

    print *, 'This program will determine the best fit parameters for the semi-empirical mass formula.'
    print *, 'The program also will write theoretical values and its uncertainty to an output file.'
    print *, 'In another file the program will write the positions of stable isotopes and drip-line positions.'
    print *, 'To accomplish this the program requires a function containing measured binding energies.'

    print *, 'please provide the file name with the experimental data'
    read(*, '(a)') filename

    ! when trying to open a file provided by the user it's good practice to
    ! check if the file exists in the current directory
    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
        ! Open the file and read the data. (don't forget to close the file)
        open(newunit = file_unit, file = filename, status = 'old') !action = 'read'

        !  * read the number of data from the first line in the file
        read(file_unit, *) n_points

        ! * Allocate Memory
        allocate(n_protons(1:n_points))
        allocate(n_neutrons(1:n_points))
        allocate(exp_values(1:n_points))
        allocate(uncertainties(1:n_points))

        !  * skip the two following lines (there's no data in those!)
        ! may need to either specify number of lines or detect it somehow
        do 30 skip_lines = 2, 3
            read(file_unit, *)
        30 continue

        ! initialize i
        i = 1

        !  * read the data and store it in the arrays you just allocated
        do i_do = 4, (n_points + 3)
            ! assign values from the data needed
            read(file_unit, *, iostat=stat) Nuc, A, N_exp, n_neutrons_r, n_protons_r, exp_values_r, EA, uncertainties_r
            if (stat == iostat_end) exit
            ! assigning values to our arrays
            n_protons(i) = n_protons_r
            n_neutrons(i) = n_neutrons_r
            exp_values(i) = exp_values_r
            uncertainties(i) = uncertainties_r
            !increment i
            i = i + 1
        enddo

        close(file_unit)
        ! There's additional kind of data in the file, make sure you only read
        ! the data you need. You need the columns labeled as N, Z, E and dE
    else
        ! Give a message saying that the file can not be found and stop the program
        print *, 'The file could not be found, PROGRAM STOP'
        stop
    endif
end subroutine read_exp_data

!-----------------------------------------------------------------------
!! Subroutine: write_predictions
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine calculates the theoretical binding energies for each
!! isotope. The subroutine will also determine the error associated with
!! each experimental binding energy. Then the a file named results.dat is
!! written, which contains the number of protons and neutrons, as well as
!! the experimental values, experimental uncertainties, theoretical values
!! and theoretical error.
!!----------------------------------------------------------------------
!! Input:
!!
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! covariance       real        Array containing the elements of the covariance matrix
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp) :: theor_val, theor_error
    integer :: file_unit, i, num_data
    character(len=*), parameter :: file_name = 'results.dat'

    ! opening file
    open(newunit=file_unit,file=file_name)

    ! writing header
    write(file_unit,*) ' Protons ', ' Neutrons ', ' Experimental Binding Energy ', &
     ' Experimental Error ', ' Theoretical Binding Energy ', ' Theoretical Error '
    num_data = size(n_protons)

    ! writing values
    do i =1, size(n_protons)
        theor_val = semi_empirical_mass(c_parameters, n_protons(i), n_neutrons(i))
        theor_error = semi_empirical_error(covariance, n_protons(i), n_neutrons(i))
        write(file_unit, '(4x, i4, 6x, i4, 10x, e13.6, 14x, e13.6, 18x, e13.6, 22x, e13.6)')&
        n_protons(i), n_neutrons(i), exp_values(i), uncertainties(i), theor_val, theor_error
    enddo
    close(file_unit)
    print *, 'theoretical binding energies were written in ', file_name
end subroutine write_predictions

!-----------------------------------------------------------------------
!! Subroutine: write_advanced
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine determines the drip-line and valley of stability positions.
!! To accomplish this the program calls the valley_of_stability and
!! theoretical_driplines subroutines and then writes the results to 
!! a file named 'advanced_results.dat'.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------

subroutine write_advanced(c_parameters)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    character(len=*), parameter :: file_name = 'advanced_results.dat'
    real(dp) :: theor_valley, theor_drip
    integer:: file_unit, i
    integer, allocatable :: N_STABLE(:), N_DRIP(:), protons(:)

    ! open file
    open(newunit=file_unit,file=file_name)

    !write header
    write(file_unit,*) ' Protons ', ' Valley of Stability Position ', ' Drip-Line Positions '

    allocate(N_STABLE(1:118))
    allocate(N_DRIP(1:118))
 
    ! calling advanced subroutines
    call valley_of_stability(c_parameters, protons, N_STABLE)
    call theoretical_driplines(c_parameters, N_DRIP)

    !writing values
    do i = 1, 118 ! this can't be right
        write(file_unit, '(6x, i4, 12x, i4, 18x, i4)')&
        protons(i), N_STABLE(i), N_DRIP(i)
    enddo
    close(file_unit)
    print *, 'The Stable Isotope and Drip-line Positions were written in ', file_name

end subroutine write_advanced
    
end module read_write