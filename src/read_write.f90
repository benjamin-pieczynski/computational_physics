!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! Benjamin Pieczynski
!!
!! Contains subroutines and functions related to reading input from the
!! user and  writing output into a text file. Responsible for writing
!! output files for each different potential.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_probability_density
!! print_energies
!! write_woods_saxon_energies
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use qm_solver, only: solve_woods_saxon
implicit none

private
public :: read_input, write_probability_density, print_energies, write_woods_saxon_energies

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Displays a message describing what the program does and the expected
!! input. After that it uses the `read_real` and `read_integer`
!! functions to assign values to the different parameters.
!!----------------------------------------------------------------------
!! Output:
!!
!! n_points     integer     number of grid points the discretized wave function
!! length       real        length of the box
!! radius       real        radius of the Woods-Saxon potential
!-----------------------------------------------------------------------
subroutine read_input(n_points, length, radius)
    implicit none
    integer, intent(out) :: n_points, radius
    real(dp), intent(out) :: length

    print *, 'BEGIN PROGRAM'
    print *, 'This program will solve the Schrodinger Equation for 3 different potentials!'
    print *, '(1) The Infinite Well, (2) The Harmonic Oscillator, and (3) Woods-Saxon.'
    print *, 'The program will provide a numerical solution for all three...'
    print *, 'It will also provide an analytic solution for (1) and (2).'
    print *, 'You will be prompted for three different values.'
    print *, 'Please enter a positive non-zero number for each value when prompted.'
    print *, ' '
    print *, 'USER INPUTS'

    ! double precision / floats
    length = read_real('box length: L')

    ! integer precision
    n_points = read_integer('the number of grid points for the discretized wave function: N')

    ! double precision / floats
    radius = read_real('Woods Saxon potential radius: r_max')

    print *, ' '
    print *, 'RUNNING PROGRAM'

end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function reads the user input for different variables to be used
!! by the program. This is for obtaining inputs of positive real numbers.
!! If the user has entered an incorrect input the function will notify 
!! the user.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! ret_val    real        A positive non negative number given by the user
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(ret_val)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a non-zero positive value for the '//trim(name)//':'

    ! Checking for a positive non-zero number
    do
        read(*,'(a)', iostat=ierror) string
        if (string /= '') then
            read(string, *, iostat=ierror) ret_val  
            ! real values will return ierror = 0, we need to check if those values are positive.
            if ( ierror == 0 ) then
                if (ret_val > 0) exit
                print *, "'"//trim(string)//"'"//' cannot be a negative or zero, please provide a real positive number'
            else 
                print *, "'"//trim(string)//"'"//' is not a number, please provide a positive real number'
            endif
        else
            print *, 'that was an empty input, please provide a non-zero positive number'
        endif
    enddo
end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function reads the user input for integers in order to determine
!! if the input is a real positive number. If the user has entered an 
!! incorrect input the function will notify the user. This check is for 
!! integers.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being requested.
!!----------------------------------------------------------------------
!! Output:
!!
!! ret_val        integer     A positive non negative number given by the user.
!-----------------------------------------------------------------------
integer function read_integer(name) result(ret_val)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a non-zero positive value for the '//trim(name)//':'

    ! Checking for a positive non-zero number
    do
        read(*,'(a)', iostat=ierror) string
        print *, ierror
        if (string /= '') then
            read(string, *, iostat=ierror) ret_val  
            if ( ierror == 0 ) then
                if (ret_val > 0) exit
                print *, "'"//trim(string)//"'"//' cannot be a negative or zero, please provide a real positive number'
            else 
                print *, "'"//trim(string)//"'"//' is not an integer, please provide a positive real number'
            endif
        else
            print *, 'that was an empty input, please provide a non-zero positive number'
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!! Subroutine: write_probability_density
!-----------------------------------------------------------------------
!! Benjamin Pieczynski
!!
!! This subroutine writes a file with the probability densities for each
!! different potential. The format is x, ground state, 1st excited, and
!! 2nd excited.
!!----------------------------------------------------------------------
!! input:
!!
!! file_name        character    a string with the output filename
!! x_vector         real        an array with the x positions
!! wave_functions   real        array of solved eigenvectors
!-----------------------------------------------------------------------
!! Output:
!!
!! file with name file_name
!-----------------------------------------------------------------------

subroutine write_probability_density(file_name, x_vector, wave_functions)
    implicit none
    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: x_vector(:), wave_functions(:,:)
    integer :: unit, i, N

    N = size(x_vector)

    open(newunit=unit,file=trim(file_name))
    write(unit, *) '   x    ', '         ground state        ', '        1st excited      ', '      2nd excited   '
    do i = 1, N
        write(unit, *) x_vector(i), wave_functions(i,1), wave_functions(i,2), wave_functions(i,3)
    end do
    close(unit)
    
end subroutine write_probability_density

!-----------------------------------------------------------------------
!! Subroutine: print_energies
!-----------------------------------------------------------------------
!! Benjamin Pieczynski
!!
!! This subroutine prints a statement that indicates the 3 lowest
!! energy states for a specified potential energy. The results for the
!! analytic and numerical solutions 
!!----------------------------------------------------------------------
!! input:
!!
!! name        character    a string with the name of the potential
!! numerical   real        an array with the 3 lowest numerical solutions
!! analytic    real        an array with the 3 lowest analytic solutions
!-----------------------------------------------------------------------
!! Output:
!!
!! print statements to terminal
!-----------------------------------------------------------------------

subroutine print_energies(name, numerical, analytic)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: numerical(:), analytic(:)
    integer :: i

    if (size(numerical) /= size(analytic)) then
        print*, "arrays size don't match in print_energies"
        stop
    endif

    print*, 'Comparing numerical and anlytic solutions in'
    print*, trim(name)
    print*, " "
    print'(a9,2a15)', 'number', 'numerical', 'analytic'
    do i = 1,3
        print *, i, numerical(i), analytic(i)
    end do
end subroutine print_energies

!-----------------------------------------------------------------------
!! Subroutine: write_woods_saxon_energies
!-----------------------------------------------------------------------
!! Benjamin Pieczynski
!!
!! This subroutine writes the energies for the Woods-Saxon potential from
!! r = 2 to r = 10. The subroutine takes the 3 lowest values for each 
!! radius of the woods saxon potential and then it writes them to an 
!! output file.
!!----------------------------------------------------------------------
!! input:
!!
!! file_name   character    a string with the name of the output file
!! n_points    integer      number of points
!! length      integer      length of sample box
!! r_min       real         minimum radius
!! r_max       real         maximum radius
!! x_vector    real         array with x postions
!-----------------------------------------------------------------------
!! Output:
!!
!! output file of name file_name
!-----------------------------------------------------------------------
subroutine write_woods_saxon_energies(file_name, N, L, r_min, r_max, x_vec)
    implicit none
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: N
    real(dp), intent(in) :: L
    real(dp), intent(in) :: x_vec(:)
    real(dp), allocatable :: energies(:), wave_functions(:,:)
    integer :: unit, i, j, r_min, r_max, R

    allocate(wave_functions(1:N,1:N))

    ! allocate
    allocate(energies(1:3))

    ! obtain 3 energies for each R and then write
    ! write to file
    open(newunit=unit,file=trim(file_name))
    write(unit, *) '     radius   ', '     1st lowest E   ', '      2nd lowest E   ', '     3rd lowest E'
    do R = r_min, r_max
        call solve_woods_saxon(N, L, R, energies, wave_functions, x_vec)
        write(unit,*) R, energies(1), energies(2), energies(3)
    end do

    close(unit)
    
end subroutine write_woods_saxon_energies


end module read_write
