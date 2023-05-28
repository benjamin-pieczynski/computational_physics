!-----------------------------------------------------------------------
!Module: mechanics
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module contains the ode functions that will be solved with the
!! 4th order Runge Kutta solver. This specifically contains the 8 coupled
!! first order ODEs that need to be solved. The module also contains a 
!! subroutine that will solve for the total energy in the system. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_energy
!! calculat_iomega
!!----------------------------------------------------------------------
!! Included functions:
!!
!! planets_ode
!-----------------------------------------------------------------------
module mechanics
use types
implicit none
private
public :: planets_ode, calculate_energy, calculate_iomega

real(dp), parameter :: G = 1.0_dp

contains


!-----------------------------------------------------------------------
!! function: planets_ode
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function generates the 8 coupled first order ordinary differential
!! equations that are needed for the runge kutta solver. The output array
!! f -> [d/dt x1, d/dt y1, d/dt x2, d/dt y2, d/dt vx1, d/dt vy1, d/dt vx2, 
!! d/dt vy2]. For simplicity G = 1.
!!----------------------------------------------------------------------
!! Input:
!!
!! r        real        array with position and velocities
!! t        real        time
!! work     real        array containing the mass for each object
!!----------------------------------------------------------------------
!! Output:
!!
!! f        real        array containing the 8 ODEs
!-----------------------------------------------------------------------
function planets_ode(r, t, work) result(f)
    implicit none
    real(dp), intent(in) :: r(:), t, work(:)
    real(dp), allocatable :: f(:)
    real(dp) :: r_1, r_2, r_12, M_p, m1, m2

    ! r -> x1, y1, x2, y2, vx1, vy1, vx2, vy2
    ! work -> M_p, m1, m2
    ! f -> [d/dt x1, d/dt y1, d/dt x2, d/dt y2, d/dt vx1, d/dt vy1, d/dt vx2, d/dt vy2]
    allocate(f(1:8))

    ! initializing masses
    M_p = work(1)
    m1 = work(2)
    m2 = work(3)

    ! initializing rs
    r_1 = sqrt(r(1)**2 + r(2)**2)
    r_2 = sqrt(r(3)**2 + r(4)**2)
    r_12 = sqrt((r(1) - r(3))**2 + (r(2) - r(4))**2)

    ! function to be solved in ode
    f(1) = r(5)
    f(2) = r(6)
    f(3) = r(7)
    f(4) = r(8)
    f(5) = (-1.0_dp * G * M_p * r(1) / r_1**3) - (G * m2 / r_12**3) * (r(1) - r(3))
    f(6) = (-1.0_dp * G * M_p * r(2) / r_1**3) - (G * m2 / r_12**3) * (r(2) - r(4))
    f(7) = (-1.0_dp * G * M_p * r(3) / r_2**3) - (G * m1 / r_12**3) * (r(3) - r(1))
    f(8) = (-1.0_dp * G * M_p * r(4) / r_2**3) - (G * m1 / r_12**3) * (r(4) - r(2))

end function planets_ode

!-----------------------------------------------------------------------
!! Subroutine: calculate_energy
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine calculates the total energy of the system for a given
!! time/step. This is based on conservation of energy and solved using the 
!! solutions to the ODEs from planets_ODE. For simplicity G = 1.
!!----------------------------------------------------------------------
!! Input:
!!
!! r        real        array with position and velocities
!! n        integer     number of iterations / step size
!! work     real        array containing the mass for each object
!!----------------------------------------------------------------------
!! Output:
!!
!! energy   real        total energy in the system for a given time
!-----------------------------------------------------------------------
subroutine calculate_energy(r, n, work, energy)
    implicit none
    real(dp), intent(in) :: r(:,:), work(:)
    real(dp), intent(out), allocatable :: energy(:)
    real(dp) :: r_1, r_2, r_12, M_p, m1, m2
    integer, intent(in) :: n
    integer :: i

    allocate(energy(1:n))

    ! initializing masses
    M_p = work(1)
    m1 = work(2)
    m2 = work(3)

    ! r -> x1, y1, x2, y2, vx1, vy1, vx2, vy2
    ! work -> M_p, m1, m2
    ! f -> [d/dt x1, d/dt y1, d/dt x2, d/dt y2, d/dt vx1, d/dt vy1, d/dt vx2, d/dt vy2]

    ! do loop
    do i = 1, n

        ! initializing rs
        r_1 = sqrt(r(1,i)**2 + r(2,i)**2)
        r_2 = sqrt(r(3,i)**2 + r(4,i)**2)
        r_12 = sqrt((r(1,i) - r(3,i))**2 + (r(2,i) - r(4,i))**2)

        ! energy calculation
        energy(i)  = 0.5 * (m1*(r(5,i)**2 + r(6,i)**2)) + 0.5 * (m2*(r(7,i)**2 + r(8,i)**2)) - &
        (G * M_p * m1 / r_1) - (G * M_p * m2 / r_2) - (G * m1 * m2 / r_12)
    end do

end subroutine calculate_energy

!-----------------------------------------------------------------------
!! Subroutine: calculate_iomega
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine calculates the total angular of the system for a given
!! time/step. Angular momentum for this system is represented by:
!! m1*(x1*vy1 - y1*vx1)+ m2*(x2*vy2 - y2*vx2)
!!----------------------------------------------------------------------
!! Input:
!!
!! r        real        array with position and velocities
!! n        integer     number of iterations / step size
!! work     real        array containing the mass for each object
!!----------------------------------------------------------------------
!! Output:
!!
!! ang_momentum   real        total energy in the system for a given time
!-----------------------------------------------------------------------
subroutine calculate_iomega(r, n, work, ang_momentum)
    implicit none
    real(dp), intent(in) :: r(:,:), work(:)
    real(dp), intent(out), allocatable :: ang_momentum(:)
    real(dp) :: r_1, r_2, r_12, M_p, m1, m2
    integer, intent(in) :: n
    integer :: i

    
    ! initializing masses
    M_p = work(1)
    m1 = work(2)
    m2 = work(3)

    allocate(ang_momentum(1:n))

    ! do loop
    do i = 1, n

        ! r -> x1, y1, x2, y2, vx1, vy1, vx2, vy2
        ! work -> M_p, m1, m2

        ang_momentum(i) = m1*(r(1,i)*r(6,i) - r(2,i)*r(5,i)) + m2*(r(3,i)*r(8,i) - r(4,i)*r(7,i))

    end do

end subroutine calculate_iomega

   
end module mechanics