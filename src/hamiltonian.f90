!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module computes the hamiltonian for each different potential.
!! The hamiltonian includes the construction of the kinetic energy and
!! the potential energy for the infinite well, harmonic oscillator, and
!! woods-saxon. These are constructed to be used as diagonals and off
!! diagonals for matrices in eigen_solver.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! construct_diagonals
!! potential_diagonal
!! woods_diagonal
!! 
!!----------------------------------------------------------------------
module hamiltonian
use types
implicit none

! hardcoding parameters
! v0 - potenial depth MeV
real(dp), parameter :: v0 = 50.0_dp
! h_bar MeV fm
real(dp), parameter :: h_bar = 197.3_dp
! mass - MeV
real(dp), parameter :: mass = 939.0_dp
! a - fm
real(dp), parameter :: a = 0.2_dp

private
public construct_diagonals, potential_diagonal, woods_diagonal

contains

!-----------------------------------------------------------------------
!! Subroutine: construct_diagonals
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine constructs the Kinetic Energy diagonal and off diagonal.
!! The diagonal follows the equation h_bar**2 / mass / dx**2. The
!! off-diagonal follows the equation -0.5_dp * h_bar**2 / mass / dx**2.
!!----------------------------------------------------------------------
!! Input:
!!
!!   N           integer                     number of points
!!   dx           real                        step size
!-----------------------------------------------------------------------
!! Output:
!!
!!  KE_diagonal         real            array describing the diagonal points
!!  KE_offdiagonal      real            array describing the off-diagonal points
!-----------------------------------------------------------------------

subroutine construct_diagonals(N, dx, KE_diagonal, KE_offdiagonal)
    implicit none
    real(dp), intent(in) :: dx
    integer, intent(in) :: N
    real(dp), intent(out), dimension(1: N) :: KE_diagonal
    real(dp), intent(out), dimension(1 : N-1) :: KE_offdiagonal

    ! constructing diagonals
    KE_diagonal = h_bar**2 / mass / dx**2
    KE_offdiagonal = -0.5_dp * h_bar**2 / mass / dx**2

end subroutine construct_diagonals

!-----------------------------------------------------------------------
!! Subroutine: potential_diagonal
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine constructs the Potential Energy diagonal for the 
!! Harmonic oscillator. This is described by the equation:
!! v_diagonal = (h_bar**2 / (2*mass)) * x_vec**2.
!!----------------------------------------------------------------------
!! Input:
!!  
!!  N                integer             number of points
!!  x_vec            real                array containing x position
!-----------------------------------------------------------------------
!! Output:
!!
!!  v_diagonal       reak               array containing the potential energies
!-----------------------------------------------------------------------

subroutine potential_diagonal(N, x_vec, v_diagonal)
    implicit none
    real(dp), intent(in) :: x_vec(:)
    integer, intent(in) :: N
    real(dp), intent(out), dimension(1: N) :: v_diagonal

    ! constructing diagonals
    v_diagonal = (h_bar**2 / (2*mass)) * x_vec**2

end subroutine potential_diagonal

!-----------------------------------------------------------------------
!! Subroutine: woods_diagonal
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine constructs the Potential Energy diagonal for the 
!! Woods-Saxon potential. To find the woods-saxon potential at each
!! point the subroutine solves the following equation:
!! V = -v0 / (1 + exp((abs(x) - R) / a))
!!----------------------------------------------------------------------
!! Input:
!!  
!!  N               integer                 number of points
!!  R               real                        radius
!!  x_vec           real                    array containing x position
!-----------------------------------------------------------------------
!! Output:
!!
!! ws_diagonal      real                    array of potentials for woods-saxon
!-----------------------------------------------------------------------
subroutine woods_diagonal(N, radius, x_vec, ws_diagonal)
    implicit none
    real(dp), intent(in) :: x_vec(:)
    integer, intent(in) :: radius
    integer, intent(in) :: N 
    real(dp), intent(out) :: ws_diagonal(1:N)

    ws_diagonal = -v0 / (1 + exp((abs(x_vec) - radius) / a))

end subroutine woods_diagonal

end module hamiltonian