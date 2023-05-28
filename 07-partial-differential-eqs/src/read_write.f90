!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module is responsible for reading the namelist, and writing
!! the values determined by the program.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_time_evolution
!! write_expectation_values
!!----------------------------------------------------------------------
module read_write
use types

implicit none

private
public :: read_input, write_time_evolution, write_expectation_values

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine checks to see if the user has placed a namelist file
!! within the command arguments. If the user has not, then the program
!! chooses to use default values
!!----------------------------------------------------------------------
!! Input:
!! None
!!----------------------------------------------------------------------
!! Output:
!!
!! length       real        the length of the well
!! n_points     integer       number of points
!! n_steps      integer       number of steps
!! delta_t      real          time steps
!! width        real          initial width
!! center       real          center location
!! k_oscillator real          K oscillator potential
!!----------------------------------------------------------------------
subroutine read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator &
    , time_file, density_file)
    implicit none
    real(dp), intent(out) :: length, delta_t, width, center, k_oscillator
    integer, intent(out) :: n_points, n_steps
    character(len=*) :: time_file, density_file 
    integer :: unit1, ios, n_arguments
    character(len=1024) :: namelist_file
    logical :: file_exists

    namelist /integration/ length, n_points, n_steps, delta_t
    namelist /wave_function/ width, center
    namelist /oscillator/ k_oscillator
    namelist /output/ time_file, density_file

    length = 5._dp
    n_points = 100
    n_steps = 100
    delta_t = 0.05_dp
    width = 0.5_dp
    center = 0.2_dp
    k_oscillator = 1._dp
    time_file = 'time_results.dat'
    density_file = 'density_results.dat'

    n_arguments = command_argument_count()
    if (n_arguments == 1) then
        ! Get namelist file name
        call get_command_argument(1, namelist_file)
        ! Check if the file exists
        inquire(file = trim(namelist_file), exist=file_exists)
        if (file_exists) then
            open(newunit=unit1, file=namelist_file)
            read(unit1, nml=integration, iostat=ios)
            if(ios /= 0) then
                print*, "Error reading integration namelist"
                stop
            endif
            read(unit1, nml=wave_function, iostat=ios)
            if(ios /= 0) then
                print*, "Error reading wave_function namelist"
                stop
            endif
            read(unit1, nml=oscillator, iostat=ios)
            if(ios /= 0) then
                print*, "Error reading oscillator namelist"
                stop
            endif
            read(unit1, nml=output, iostat=ios)
            if(ios /= 0) then
                print*, "Error reading output namelist"
                stop
            endif
        ! check to see if namelist file doesn't exist
        else
            print*, namelist_file, 'not found'
            stop
        endif
    ! Error: more than one argument, stop the program
    elseif (n_arguments/=0) then
        print*, 'Not enought arguments. Program requires either 0 or 1 arguments'
        stop
    endif

end subroutine read_input


!-----------------------------------------------------------------------
!! Subroutine: write_time_evolution
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine writes to the density_file file the probability 
!! density at different times. The first LINE should contain the sample 
!! points along the x axis. The successive lines contains the 
!! probability density at different time steps.
!!----------------------------------------------------------------------
!! Input:
!!
!! N                integer     number of points
!! x_vector         real        vector containing position of lattice
!! time_wave_func   real        matrix of time_wave function
!! output_file      character   name of the output file
!!----------------------------------------------------------------------
subroutine write_time_evolution(nstep, x_vector, time_wave_func, output_file)

    !This subroutine should write to the density_file file the probability 
    !density at different times. The first LINE should contain the sample 
    !points along the x axis.

    !The successive lines should contain the probability density at  
    !different time steps.
    implicit none
    real(dp), intent(in) :: x_vector(:), time_wave_func(:,:)
    character(len=*), intent(in) :: output_file
    integer, intent(in) :: nstep
    integer :: i
    integer :: file_unit
    
    ! Open the output file
    open(newunit=file_unit, file=output_file)
    
    write(file_unit, *) x_vector(:)
    do i = 1, nstep + 1
        write(file_unit, *) time_wave_func(:,i)
    end do
    close(file_unit)
    print *, shape(time_wave_func)
    print *, size(x_vector)
    
end subroutine write_time_evolution

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine writes to the time_file file the expectation values 
!! as a function time. The first COLUMN contains the times at which the 
!! wave function was calculated. The successive columns contain the 
!! expectation values (normalization, position, width) at the respective 
!! times.
!!----------------------------------------------------------------------
!! Input:
!!
!! N                integer     number of points
!! time_matrix         real     matrix that contains time values
!! norm                real     vector containing normalization
!! position            real     vector containing positions
!! width               real     vector containing the width
!! output_file      character   name of the output file
!!----------------------------------------------------------------------
subroutine write_expectation_values(nstep, dt, norm, position, width, output_file)
    implicit none
    real(dp), intent(in) :: norm(:), position(:), width(:), dt
    character(len=*), intent(in) :: output_file
    integer, intent(in) :: nstep
    integer :: i, j
    integer :: file_unit
    real(dp) :: time

    ! Open the output file
    open(newunit=file_unit, file=output_file)
    ! Header
    write(file_unit, 40) 't', 'norm', 'pos', 'width'
    40 format(15X, A1, 15X, A8, 15X, A8, 15X, A8)

    ! DO loop: compute phi(x,y) and write to file
    do i=1, nstep+1
        time = (i * dt) - dt
        write(file_unit, *) time, norm(i), position(i), width(i)
    !   print *, time_matrix(:,i)
    enddo
    ! write(file_unit, 50) time_matrix(1,:), norm(:), position(:), width(:)
    !*50 format(5X, F7.5, 5X, F7.5, 5X, F7.5, 5X, F7.5)
    
end subroutine write_expectation_values

end module read_write
