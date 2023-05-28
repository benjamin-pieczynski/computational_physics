!-----------------------------------------------------------------------
!Module: qm_solver
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module solves for the eigenvectors and eigenvalues for each different
!! potential using the schrodinger equation. This includes the infinite
!! well(particle in a box), Harmonic oscillator, and Woods-Saxon. This 
!! module includes subroutines that construct the diagonals for Kinetic
!! and potential energy, and then passes them to eigen_solver in order
!! to determine the eigenvalues and eigenvectors. The subroutine also
!! has an analytic solution hardcoded for comparison to the numerical
!! solutions.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_infinite_well
!! sample_box
!! analytic_in
!!----------------------------------------------------------------------
module qm_solver
use types
use hamiltonian, only: construct_diagonals, potential_diagonal, woods_diagonal
use eigen_solver, only : solve_eigenproblem
implicit none

private
public solve_infinite_well, sample_box, analytic_infinite_well, solve_harmonic_oscillator, &
    analytic_harmonic_oscillator, solve_woods_saxon, normalize_wf

real(dp), parameter :: mass = 939.0_dp
real(dp), parameter :: h_bar = 197.3_dp

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine determine the x array points for the wavefunction box.
!! This is determined by x(i) = -L + dx*(i-1), where dx is the step size
!! dx = (2 * L) / (N -L).
!!----------------------------------------------------------------------
!! Input:
!!  length   real    half of the box length
!!  n_points    integer    # of points
!-----------------------------------------------------------------------
!! Output:
!!  x_arr   real array of points for the wavefunction box
!-----------------------------------------------------------------------
subroutine sample_box(L, N, x_arr)
    implicit none
    real(dp), intent(in) :: L
    integer, intent(in) :: N
    real(dp) :: step_size
    integer :: i
    real(dp), allocatable, intent(out):: x_arr(:)

    allocate(x_arr(1:N))

    ! initialize i
    i = 0

    ! finding step size
    step_size = (2 * L) / (N - 1)

    ! do loop to find array for x
    do i = 1, N - 1
        x_arr(i) = -L + step_size*(i-1)
    end do
    ! must be L
    x_arr(N) = L

end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: solve_infinite_well
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine provides a numerical solution for the infinite well/
!! particle in a box problem. This function has potential energy set
!! to zero, so it only constructs a diagonal for the kinetic energy.
!! It then uses an eigensolver to determine the 
!!----------------------------------------------------------------------
!! Input:
!!
!!  N           integer     number of points
!!  L           real        length of the box
!!  x_vec       real        array containing x position
!-----------------------------------------------------------------------
!! Output:
!!
!!  energies    real        numerical solution with 3 lowest energies
!!  low_wfunc   real        numerical solution with 3 lowest wavefunctions
!-----------------------------------------------------------------------
subroutine solve_infinite_well(N, L, energies, low_wvfunc, x_vec)
    implicit none
    integer, intent(in) :: N ! # of points
    real(dp), intent(in) :: L ! length of the box
    real(dp), allocatable, intent(out) :: low_wvfunc(:,:)
    real(dp), intent(out) :: energies(:)
    real(dp) :: KE_diagonal(1:N)
    real(dp) :: KE_offdiagonal(1:N-1)
    real(dp) :: eigenvalues(1:N)
    real(dp) :: eigenvector(1:N,1:N)
    real(dp) :: dx
    real(dp), intent(in) :: x_vec(:)
    integer :: i

    allocate(low_wvfunc(1:N, 1:N))

    ! dx is just the difference between 2 x points
    dx = x_vec(2) - x_vec(1)
    call construct_diagonals(N, dx, KE_diagonal, KE_offdiagonal)
    call solve_eigenproblem(KE_diagonal, KE_offdiagonal, eigenvalues, eigenvector)

    ! normalize eigen wavefunction
    call normalize_wf(eigenvector, dx, N)

    ! 3 lowest energies
    do i = 1, 3
        energies(i) = eigenvalues(i)
        low_wvfunc(:,i) = (eigenvector(:,i))**2
    end do

end subroutine solve_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: analytic_infinite_well
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine uses an analytic equation to solve for the infinite
!! well / particle in a box problem. It returns the analytic solution
!! for the 3 lowest energies.
!!----------------------------------------------------------------------
!! Input:
!!
!! L                      real    length of the box
!-----------------------------------------------------------------------
!! Output:
!!
!! analytic_energies      real    vector with 3 lowest energies
!-----------------------------------------------------------------------
subroutine analytic_infinite_well(L, analytic_energies)
    implicit none
    real(dp), intent(out) :: analytic_energies(:)
    real(dp), intent(in) :: L 
    integer :: i

    do i = 1, 3
        analytic_energies(i) = i**2 * h_bar**2 * pi**2 / 8 / mass / L**2
    end do
end subroutine analytic_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: solve_harmonic_oscillator
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine provides the numerical solution for the case of the
!! harmonic oscillator. It computes the kinetic energy diagonal and
!! offdiagonal as well as the potential energy diagonal through different
!! subroutines. Then it passes the diagonals and off diagonals to the
!! eigen solvver and then outputs the 3 lowest energies and wavefunctions.
!!----------------------------------------------------------------------
!! Input:
!!
!!  N           integer                  number of points
!!  L           real                     length of the box
!!  x_vec       real                     array of x positions
!-----------------------------------------------------------------------
!! Output:
!!
!!  energies        real                array with 3 lowest energies
!!  low_wvfunc      real                array with 3 lowest wavefunctions
!-----------------------------------------------------------------------
subroutine solve_harmonic_oscillator(N, L, energies, low_wvfunc, x_vec)
    implicit none
    integer, intent(in) :: N ! # of points
    real(dp), intent(in) :: L ! length of the box
    real(dp), intent(out) :: low_wvfunc(:,:), energies(:)
    real(dp) :: KE_diagonal(1:N)
    real(dp) :: KE_offdiagonal(1:N-1)
    real(dp) :: v_diagonal(1:N)
    real(dp) :: eigenvalues(1:N)
    real(dp) :: eigenvector(1:N,1:N)
    real(dp), intent(in) :: x_vec(:)
    real(dp) :: dx
    integer :: i

    dx = x_vec(2) - x_vec(1)
    call construct_diagonals(N, dx, KE_diagonal, KE_offdiagonal)
    call potential_diagonal(N, x_vec, v_diagonal)
    call solve_eigenproblem(KE_diagonal + v_diagonal, KE_offdiagonal, eigenvalues, eigenvector)

    ! normalize eigen wavefunction
    call normalize_wf(eigenvector, dx, N)

    ! 3 lowest energies
    do i = 1, 3
        energies(i) = eigenvalues(i)
        low_wvfunc(:,i) = (eigenvector(:,i))**2
    end do
end subroutine solve_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: analytic_harmonic_oscillator
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine provides the analytic solution for the harmonic oscillator.
!! This is described as E = (i - 0.5) * h_bar**2 / mass
!! The subroutine outputs the 3 lowest energies for this solution.
!!----------------------------------------------------------------------
!! Input:
!!
!! none
!-----------------------------------------------------------------------
!! Output:
!!
!! harmonic_enerfy      real        array with the 3 lowest energies
!-----------------------------------------------------------------------
subroutine analytic_harmonic_oscillator(harmonic_energy)
    implicit none
    real(dp), intent(out) :: harmonic_energy(1:3)
    integer :: i

    do i = 1,3
        harmonic_energy(i) = (i - 0.5) * h_bar**2 / mass
    end do

end subroutine analytic_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: solve_woods_saxon
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine provides the numerical solution for the case of the
!! woods-saxons potential. It computes the kinetic energy diagonal and
!! offdiagonal as well as the potential energy diagonal through different
!! subroutines. Then it passes the diagonals and off diagonals to the
!! eigen solvver and then outputs the 3 lowest energies and wavefunctions.
!!----------------------------------------------------------------------
!! Input:
!!
!!  N           integer                  number of points
!!  L           real                     length of the box
!!  x_vec       real                     array of x positions
!!  R           integer                     woods-saxon radius
!-----------------------------------------------------------------------
!! Output:
!!
!!  energies        real                array with 3 lowest energies
!!  low_wvfunc      real                array with 3 lowest wavefunctions
!-----------------------------------------------------------------------
subroutine solve_woods_saxon(N, L, R, energies, low_wvfunc, x_vec)
    implicit none
    integer, intent(in) :: R, N ! # of points
    real(dp), intent(in) :: L
    real(dp), intent(out) :: low_wvfunc(:,:), energies(:)
    real(dp) :: KE_diagonal(1:N)
    real(dp) :: KE_offdiagonal(1:N-1)
    real(dp) :: v_diagonal(1:N)
    real(dp) :: eigenvalues(1:N)
    real(dp) :: eigenvector(1:N,1:N)
    real(dp), intent(in) :: x_vec(:)
    real(dp) :: dx
    integer :: i
    
    dx = x_vec(2) - x_vec(1)
    call construct_diagonals(N, dx, KE_diagonal, KE_offdiagonal)
    call woods_diagonal(N, R, x_vec, v_diagonal)
    call solve_eigenproblem(KE_diagonal + v_diagonal, KE_offdiagonal, eigenvalues, eigenvector)
    
    ! normalize eigen wavefunction
    call normalize_wf(eigenvector, dx, N)
    
    ! 3 lowest energies
    do i = 1, 3
        energies(i) = eigenvalues(i)
        low_wvfunc(:,i) = (eigenvector(:,i))**2
    end do

end subroutine solve_woods_saxon

!-----------------------------------------------------------------------
!! Subroutine: normalize_wf
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine normalizes incoming eigen vectors / wave functions
!! and then returns the normalized function.
!!----------------------------------------------------------------------
!! Input:
!!
!!  eigen_vector        real             eigen vector to normalize
!!  d_x                 real             step size
!!  m                   integer          number to normalize
!-----------------------------------------------------------------------
!! Output:
!!
!!  eigen_vector       real              normalized wavefunctions
!-----------------------------------------------------------------------
subroutine normalize_wf(eigen_vector, d_x, m)
    implicit none
    real(dp), intent(inout) :: eigen_vector(:,:)
    real(dp), intent(in) :: d_x
    integer, intent(in) :: m
    integer :: i

    do i = 1, m
        eigen_vector(:,i) = eigen_vector(:,i)*(1._dp/sqrt(d_x*sum(eigen_vector(:,i)**2)))
    end do

end subroutine normalize_wf

end module qm_solver