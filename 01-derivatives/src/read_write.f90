!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module is responsible for calling the 3 and 5 point functions in
!! euler_formulas.f90 and looping them until h_step > h_max. The
!! analytic solution function, analytic_f, is also called from the
!! module analytic_functions.f90. This module is also responsible for 
!! asking for a user input for x_zero and checking that the user has 
!! entered a real number.  
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input(x_zero)
!! write_derivatives(x_zero)
!-----------------------------------------------------------------------
module read_write
use types
use analytic_functions, only : second_derivative_f
use euler_formulas, only : euler_3points, euler_5points
implicit none

! The private statement restricts every function, parameter and variable
! defined in this module to be visible only by this module
private
! Then we use the public statement to only make visible to other modules 
! the few functions or subroutines that will be used by them
public read_input, write_derivatives

contains

!-----------------------------------------------------------------------
!Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine informs the user that the program calculates the
!! second derivative of the function x*sin(x) using 3 and 5 point
!! euler derivative approximations.  The subroutine then requsts that
!! the user input a real number for x_zero in order to calculate the
!! 2nd derivative.  After recieving the input the subroutine will check
!! to make sure that the user has inputted by reading a string of the
!! user's input.  If the input is empty or does not contain a real
!! number the program will notify the user and request a real number.
!! 
!!----------------------------------------------------------------------
!! Arguments:
!! x_zero  real  value at which the derivatives will be calculated
!-----------------------------------------------------------------------
subroutine read_input(x_zero)
    implicit none
    real(dp), intent(out) :: x_zero
    character(len=120) :: string
    integer :: ierror

    print *, 'This program calculates the second derivative of x*sin(x)'
    print *, 'The calculation uses both the 3 and 5 point euler derivative approximations'
    print *, 'The results for both the 3 point and 5 point methods will be available in output.txt'
    print *, 'In order to run this program please input at real number for x_zero:'


    ! We enclose the input reading inside an infinite loop that can only
    ! be exited when a correct input is given.
    !
    ! Instead of trying to read a real number we read a string containing
    ! the user's input and then make checks on that string by converting 
    ! it into a real number.
    ! 
    ! The first check is to make sure that the string is not empty 
    ! (i.e. the user simply pressed the enter key)
    ! 
    ! The second check is made by using the 'read' statement to convert
    ! the string into a number, if that is not possible iostat gives an
    ! error code different from zero.
    do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x_zero
            if (ierror == 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo
end subroutine read_input

!-----------------------------------------------------------------------
!Subroutine: write_derivatives
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine assigns values to h_step, h_increment and h_max.
!! These variables will be utilized by both the 3 and 5 point
!! formulas defined in analytic_functions.f90 to calculate the second
!! derivative of x*sin(x).  
!! 
!!----------------------------------------------------------------------
!! Arguments:
!! x_zero  real  value at which the derivatives will be calculated
!-----------------------------------------------------------------------
subroutine write_derivatives(x_zero)
    implicit none
    real(dp), intent(in) :: x_zero
    real(dp) :: h_step = 0.2_dp
    real(dp), parameter :: h_increment = 1.01_dp, h_max = 0.57_dp
    real(dp) :: d2_analytic, d2_num3, d2_num5
    character(len=*), parameter :: file_name = 'results.dat'
    integer :: unit

    d2_analytic = second_derivative_f(x_zero)

    open(newunit=unit, file=file_name)
    write(unit,'(4a28)') 'h', 'analytic', '3 point', '5 point'
    do 
        d2_num3 = euler_3points(x_zero, h_step)
        d2_num5 = euler_5points(x_zero, h_step)
        write(unit,'(4e28.16)') h_step, d2_analytic, d2_num3, d2_num5
        if(h_step > h_max) exit
        h_step = h_step * h_increment
    enddo
    close(unit)

    print *, 'The derivatives were written in the '//file_name//' file'
end subroutine write_derivatives

end module read_write
