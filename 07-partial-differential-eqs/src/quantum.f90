!-----------------------------------------------------------------------
!Module: quantum
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module is responsible for computing the expectation values and
!! for handling the time evolution matrix. This module contains the 
!! Crank-Nicolson method, where we split the hamiltonian into half
!! explicit and half implicit. We use a method to avoid using imaginary
!! values in order solve for the time dependent schrodinger equation.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! sample_box
!! construct_initial_wavefunction
!! evolve_wave_function
!! construct_time_evolution_matrix
!! expectation_values
!!----------------------------------------------------------------------

module quantum

use types
use hamiltonian, only: construct_diagonals, potential_diagonal
use linear_algebra, only: invert_matrix

implicit none
private
public sample_box, construct_initial_wavefunction, evolve_wave_function, &
construct_time_evolution_matrix, expectation_values

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
!!  L  real    half of the box length
!!  N    integer    # of points
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

    ! finding step size
    step_size = (2._dp * L) / (N - 1)

    ! do loop to find array for x
    ! might need to initialize x_arr = 0._dp
    do i = 1, N
        x_arr(i) = -L + step_size*(i-1)
    end do

end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: construct_initial_wavefunction
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
!!  center      real        center position
!!  width       real        specified width
!!  x_vec       real        array containing x position
!-----------------------------------------------------------------------
!! Output:
!!
!!  width    real        numerical solution with 3 lowest energies
!!  wvfunc   real        wavefunction half implicit, half explicit
!-----------------------------------------------------------------------
subroutine construct_initial_wavefunction(N, center, width, wvfunc, x_vec)
        implicit none
        integer, intent(in) :: N ! # of points
        real(dp), intent(in) :: center, width
        real(dp), allocatable, intent(out) :: wvfunc(:)
        real(dp), allocatable, intent(in) :: x_vec(:)
        integer :: i
    
        allocate(wvfunc(1:2*N))
    
        do i = 1, 2*N
            if (i .GT. N) then
                wvfunc(i) = 0._dp
            else
                wvfunc(i) = ((2._dp*pi*width**2._dp)**(-1._dp/4._dp)) * & 
                (exp((x_vec(i)-center)**2._dp / (-4._dp*width**2._dp)))
            end if
        end do

end subroutine construct_initial_wavefunction

!-----------------------------------------------------------------------
!! Subroutine: construct_time_evolution_matrix
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine constructs the time evolution matrix.
!!----------------------------------------------------------------------
!! Input:
!!
!!  N           integer     number of points
!!  L           real        length of box
!!  delta_t     real        time spacing
!!  x_vec       real        array containing x position
!!----------------------------------------------------------------------
!! Output:
!!
!! evolution_matrix     real    matrix containing the time evolution
!!----------------------------------------------------------------------
subroutine construct_time_evolution_matrix(N, L, k_oscillator, x_vec, delta_t, evolution_matrix)
    implicit none
    real(dp), intent(in) :: delta_t, L, k_oscillator
    real(dp), allocatable, intent(in) :: x_vec(:)
    real(dp), allocatable, intent(out) :: evolution_matrix(:,:)
    real(dp), allocatable :: KE_diagonal(:), KE_offdiagonal(:), v_diagonal(:), H_diagonal(:), &
    alpha(:,:), beta(:,:), a_m1(:,:), H_matrix(:,:), identity_matrix(:,:)
    integer, intent(in) :: N
    real(dp) :: mods
    integer :: i, j

    ! memory allocation
    allocate(evolution_matrix(1:2*N,1:2*N))
    allocate(alpha(1:2*N,1:2*N))
    allocate(beta(1:2*N,1:2*N))
    allocate(H_diagonal(1:N))
    allocate(H_matrix(1:N,1:N))
    allocate(identity_matrix(1:N,1:N))

    ! if k_osc = 0 the basic portion runs, set to a value for advanced
    call potential_diagonal(N, k_oscillator, x_vec, v_diagonal)
    mods = delta_t / (2._dp*h_bar)

    call construct_diagonals(N, L, KE_diagonal, KE_offdiagonal)

    ! hamiltonian matrix
    h_diagonal(:) = KE_diagonal(:) + v_diagonal(:)
    H_matrix(:,:) = 0._dp
    do i = 1, N
        H_matrix(i,i) = H_diagonal(i)
        if (i .LT. N) then
            H_matrix(i,i+1) = KE_offdiagonal(i)
            H_matrix(i+1,i) = KE_offdiagonal(i)
        end if
    end do

    ! diagnal is 1 for both matrices (upper left, lower right quadrants)
    ! identity matrix
    identity_matrix = 0._dp
    do i = 1, N
        identity_matrix(i,i) = 1._dp
    end do

    ! initialize our matrices to 0
    alpha(:,:) = 0._dp
    beta(:,:) = 0._dp

    ! constructing the alpha and beta matrixes
    do i = 1, N
        do j = 1, N  
            alpha(i,j) = identity_matrix(i,j)
            beta(i,j) = identity_matrix(i,j)
            alpha(i+N,j+N) = identity_matrix(i,j)
            beta(i+N,j+N) = identity_matrix(i,j)
            alpha(i+N,j) = -mods * H_matrix(i,j)
            beta(i,j+N) = -mods * H_matrix(i,j)
            beta(i+N,j) = mods * H_matrix(i,j)
            alpha(i,j+N) = mods * H_matrix(i,j)
        end do
    end do

    ! invert the alpha matrix
    call invert_matrix(alpha, a_m1)

    ! construct the time evolution matrix
    evolution_matrix = matmul(a_m1, beta)

end subroutine construct_time_evolution_matrix

!-----------------------------------------------------------------------
!! Subroutine: evolve_wave_function
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function is responsible for evolving the wave function over time.
!!----------------------------------------------------------------------
!! Input:
!!
!!  N           integer     number of points
!!  nstep       integer     number of steps
!!  L           real        length of box
!!  wave_func   real        the wave function (array)
!!  evolution_matrix real   real time evolution matrix
!!----------------------------------------------------------------------
!! Output:
!!
!! time_wave_func   real      contains the density wave function
!!----------------------------------------------------------------------
subroutine evolve_wave_function(N, nstep, wave_func, evolution_matrix, time_wave_func)
    implicit none
    real(dp), allocatable, intent(in) :: evolution_matrix(:,:), wave_func(:)
    real(dp), allocatable, intent(out) :: time_wave_func(:,:)
    !real(dp), intent(in) :: L
    real(dp), allocatable :: rho(:,:)
    integer, intent(in) :: N, nstep
    integer :: i, j
    !real(dp) :: dx

    !dx = 2*L / (N - 1)

    ! allocate memory
    allocate(time_wave_func(1:N,1:(nstep+1)))
    allocate(rho(1:(2*N),1:(nstep+1)))

    rho(:,1) = wave_func(:)
    do i = 1, nstep
        rho(:,i+1) = matmul(evolution_matrix, rho(:,i))
    end do 

    ! getting times
    do i = 1, nstep + 1 ! we need the first t=0
        do j = 1, N
            time_wave_func(j,i) = (rho(j,i)**2.0_dp) + (rho(j+N,i)**2.0_dp)
        end do
    end do
    
end subroutine evolve_wave_function

!-----------------------------------------------------------------------
!! Subroutine: expectation_values
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine determines the expectation values for the normalization,
!! the position, and the widths. This is for the time dependent 
!! Schrodinger equation.
!!----------------------------------------------------------------------
!! Input:
!!
!! time_wave_func   real      contains the density wave function
!!  N           integer     number of points
!!  n_steps      integer     number of steps
!!  wave_func   real        the wave function (array)
!!  L           real        length of box
!!  x_vec       real        array containing x position
!!----------------------------------------------------------------------
!! Output:
!!
!! norm                real     vector containing normalization
!! position            real     vector containing positions
!! sigma              real     vector containing the width
!!----------------------------------------------------------------------
subroutine expectation_values(N, L, n_steps, x_vec, time_wave_func, norm, position, sigma)
    implicit none
    real(dp), intent(in) :: time_wave_func(:,:), x_vec(:), L
    real(dp), allocatable, intent(out) :: norm(:), position(:), sigma(:)
    real(dp) :: x_2, dx, t_norm, x, t_sigma
    integer, intent(in) :: N, n_steps
    integer :: i

    ! allocate memory
    allocate(norm(1:n_steps+1))
    allocate(position(1:n_steps+1))
    allocate(sigma(1:n_steps+1))

    dx = 2._dp*L / (N - 1)*1._dp

    ! getting our expectation values
    do i = 1, n_steps + 1
        t_norm = sum(time_wave_func(:,i)*dx)
        x = sum(x_vec*time_wave_func(:,i)*dx) / t_norm
        x_2 = sum( (x_vec**2.0_dp)*time_wave_func(:,i)*dx) / t_norm
        t_sigma = sqrt(x_2 - x**2._dp)
        norm(i) = t_norm ! equal to temp val
        sigma(i) = t_sigma ! equal to temp val
        position(i) = x
    end do

end subroutine expectation_values
end module quantum