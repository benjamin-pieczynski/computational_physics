!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module reads the input for possible namelist arguments and then
!! writes the results into an output file after the program has determined
!! the solution to the 8 coupled ODEs. If a namelist argument is detected
!! the program will pull the initial values from the namelist, if not the
!! values will be taken from the program default values. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_results
!!----------------------------------------------------------------------

module read_write
use types
implicit none

private
public :: read_input, write_results

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine looks for a user namelist argument. If a namelist 
!! argument is detected the program will pull the initial values from 
!! the namelist, if not the values will be taken from the program's
!! default values.
!!----------------------------------------------------------------------
!! Output:
!!
!! work_array           real        array to store masses
!! initial_condition    real        array to store initial conditions
!! final_time           real        last time in time array
!! n_steps              integer     number of steps
!! output_file          character   name of the output file
!-----------------------------------------------------------------------
subroutine read_input(work_array, initial_condition, final_time, n_steps, output_file)
    implicit none
    real(dp), intent(out) :: work_array(1:3), initial_condition(1:8)
    real(dp), intent(out) :: final_time
    integer, intent(out) :: n_steps
    integer :: n_arguments, file_unit, ierror
    logical :: file_exists
    character(len=*) :: output_file
    real(dp) :: primary_mass, planet_mass_1, planet_mass_2
    real(dp) :: initial_pos_1(1:2), initial_pos_2(1:2)
    real(dp) :: initial_vel_1(1:2), initial_vel_2(1:2)

    namelist /masses/ primary_mass, planet_mass_1, planet_mass_2
    namelist /initial_conditions/ initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2
    namelist /solution_parameters/ final_time, n_steps
    namelist /output/ output_file

    ! Set default values
    primary_mass  = 1
    planet_mass_1 = 1
    planet_mass_2 = 1
    initial_pos_1 = [1._dp, 0._dp]
    initial_pos_2 = [0._dp, 1._dp]
    initial_vel_1 = [0._dp, 1._dp]
    initial_vel_2 = [0._dp, 0.5_dp]
    final_time = 100
    n_steps = 1000
    output_file = 'planet_default.dat'

    ! get namelist file name from command line
    n_arguments = command_argument_count() 

    ! read namelists
    if (n_arguments == 1) then
        call get_command_argument(1, output_file)
        inquire(file = trim(output_file), exist = file_exists)
        if (file_exists) then
            open(newunit=file_unit, file=output_file)
            read(file_unit, nml=masses, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading masses namelist"
                stop
            endif
            read(file_unit, nml=initial_conditions, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading initial_conditions namelist"
                stop
            endif
            read(file_unit, nml=solution_parameters, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading solution_parameters namelist"
                stop
            endif
            read(file_unit, nml=output, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading output namelist"
                stop
            endif
        else
            print*, output_file, 'not found'
            stop
        endif
    elseif (n_arguments /= 0) then
        print*, 'Incorrect number of arguments. Program takes either 0 or 1 argument only'
        print*, 'See details in README.md'
        stop
    endif

    ! initializing arrays
    work_array = [primary_mass, planet_mass_1, planet_mass_2]
    initial_condition = [initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2]

end subroutine read_input

!-----------------------------------------------------------------------
!! Subroutine: write_results
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine writes the results for the program. The program writes
!! to either the default output file or a file specified by a user
!! namelist argument. Order of result variables is t, x1, y1, x2, y2, and E.
!!----------------------------------------------------------------------
!! Input:
!!
!! t            real        time array
!! r            real        solution array
!! energy       real        total energy array
!! file_name    character   output file name
!-----------------------------------------------------------------------
subroutine write_results(t, r, energy, ang_momentum, n, file_name)
    implicit none
    real(dp), intent(in) :: t(:), r(:,:), energy(:), ang_momentum(:)
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: n
    integer :: unit, i

    
    ! writing to the file
    print *, file_name
    open(newunit=unit, file=trim(file_name))
    write(unit,*) '      t    ', '                x1        ', '           y1      ',&
    '            x2      ', '         y2      ', '           E    ', '     I      '
    do i = 1, n
        write(unit,'(7e25.8)') t(i), r(1,i), r(2,i), r(3,i), r(4,i), energy(i), ang_momentum(i)
    end do
    close(unit)
end subroutine write_results


end module read_write
