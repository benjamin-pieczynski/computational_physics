!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module computes the hamiltonian for each different potential.
!! The hamiltonian includes the construction of the kinetic energy and
!! the potential energy for the infinite well, harmonic oscillator, and
!! woods-saxon. These are constructed to be used as diagonals and off
!! diagonals for matrices in eigen_solver. This has been modified
!! for only the harmonic and infinite well
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! construct_diagonals
!! potential_diagonal
!! 
!!----------------------------------------------------------------------
module hamiltonian
use types
implicit none

! hardcoding parameters
! v0 - potenial depth MeV
!real(dp), parameter :: v0 = 50.0_dp
! h_bar MeV fm
!real(dp), parameter :: h_bar = 1._dp
! mass - MeV
!real(dp), parameter :: mass = 1._dp
! a - fm
!real(dp), parameter :: a = 0.2_dp

private
public construct_diagonals, potential_diagonal

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

subroutine construct_diagonals(N, L, KE_diagonal, KE_offdiagonal)
    implicit none
    real(dp) :: dx
    real(dp), intent(in) :: L
    integer, intent(in) :: N
    real(dp), allocatable, intent(out) :: KE_diagonal(:), KE_offdiagonal(:)

    dx = (2._dp*L) / (N - 1)*1._dp
    ! constructing diagonals
    allocate(KE_diagonal(1:N))
    allocate(KE_offdiagonal(1:N-1))

    KE_diagonal(:) = h_bar**2._dp / mass / dx**2._dp
    KE_offdiagonal(:) = -0.5_dp * h_bar**2._dp / mass / dx**2._dp

end subroutine construct_diagonals

!-----------------------------------------------------------------------
!! Subroutine: potential_diagonal
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine constructs the Potential Energy diagonal for the 
!! Harmonic oscillator. This is described by the equation:
!! v_diagonal(i) = 0.5_dp * k * x(i)**2._dp
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

subroutine potential_diagonal(N, k, x_vec, v_diagonal)
    implicit none
    real(dp), intent(in) :: x_vec(:), k
    integer, intent(in) :: N
    real(dp), allocatable, intent(out) :: v_diagonal(:)
    integer :: i

    allocate(v_diagonal(1:N))
    ! constructing diagonals

    do i = 1, N
        v_diagonal(i) = 0.5_dp * k * x_vec(i)**2._dp
    end do

end subroutine potential_diagonal

end module hamiltonian