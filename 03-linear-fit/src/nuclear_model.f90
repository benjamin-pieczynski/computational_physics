!-----------------------------------------------------------------------
!Module: nuclear_model
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module determines the theoretical values and uncertainties for
!! isotopes over different N and Z values. The module also constructs
!! the alpha and beta matricies that will be used to determine our best
!! fit patameters for the theoretical model. To accomplish this we consider
!! the semi-empirical mass formula to estimate the drip line location, and
!! to determine the values for c_volume and c_surface. This builds the best-fit
!! thoretical model that is then used to calculate the binding energy for the isotopes. 
!! The module ultimately will determine the most stable isotope for Z values between 
!! 1 and 118 to find  the posions of valley of stability, as well as the neutron 
!! drip-line positions for these Z-values(protons).
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! find_best_parameters
!! construct_alpha_beta
!! calculate_linear_terms
!! valley_of_stability
!! theoretical_driplines
!! 
!!----------------------------------------------------------------------
!! Included functions:
!!
!! volume_term
!! surface_term
!! coulomb_term
!! surface_term
!! asymmetry_term
!! pairing_term
!! initialize_delta_ZN
!! semi_empirical_mass
!! semi_empirical_error
!! energy_per_nucleon
!!
!!----------------------------------------------------------------------
module nuclear_model
use types
use linear_algebra, only : solve_linear_system
implicit none

private

public :: find_best_parameters, semi_empirical_mass, semi_empirical_error, &
valley_of_stability, theoretical_driplines
contains


!-----------------------------------------------------------------------
!! Subroutine: find_best_parameters
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine finds the best fit parameters and the uncertainty for 
!! each parameter. The parameters this subroutine is finding best-fits
!! for is c_vol, c_surf, c_asym, c_coul, and c_pair. This is done by solving
!! a linear system of equations using matricies. The subroutine constructs
!! an alpha and beta matrix(calls a subroutine), inverts the alpha matrix 
!! and then solves for the linear equation to find each parameter.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! c_parameters     real        Array containing the semi-empirical mass formula parameters
!! covariance       real        Array containing the covariance matrix of the parameters
!-----------------------------------------------------------------------
subroutine find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out), allocatable ::  c_parameters(:), covariance(:,:)

    ! If you do the extra credit you would change this to 6
    integer, parameter :: n_parameters = 5

    real(dp) :: alpha(1:n_parameters,1:n_parameters), beta(1:n_parameters)

    ! c_parameters and covariance were passed as arguments and need to be
    ! allocated. Deallocate them if they're allocated and then allocate them
    ! with the correct size
    allocate(c_parameters(1:n_parameters))
    allocate(covariance(1:n_parameters, 1:n_parameters))

    call construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)

    ! The subroutine below (defined in the linear_algebra module) should solve
    ! the matrix equation in the README and return the c_parameters and the
    ! inverse of alpha (the covariance matrix)
    call solve_linear_system(alpha,beta, c_parameters, covariance)
    ! Now just print the parameters (with it's uncertainties) to screen
    call print_best_parameters(c_parameters, covariance)
end subroutine find_best_parameters

!-----------------------------------------------------------------------
!! Subroutine: construct_alpha_beta
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine is responsible for constructing the alpha and beta
!! matrix needed to find the best parameters. The subroutine first checks
!! to make sure that the size of the matrices to be constructed (1) - is
!! a square matrix if it is the alpha matrix, and (2) that beta matrix
!! has the same number of elements as the alpha matrix has rows or columns.
!! Then both matrices are constructed within a do loop that correctly sizes
!! them based on the arrays passed into this subroutine.
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! alpha            real        Array containing the alpha matrix
!! beta             real        Array containing the beta vector
!-----------------------------------------------------------------------
subroutine construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out) :: alpha(:,:), beta(:)

    integer :: n_data, n_parameters, i, j, k, i_size, j_size
    real(dp):: linear_terms(1:size(beta))

    ! size of matrix A, B
    n_parameters = size(beta)

    ! Check if the alpha array is a square matrix
    if (size(alpha, 1) /= size(alpha, 2)) then
        print *, 'ERROR: alpha is not a square matrix'
        stop
    endif

    ! Also check that beta has the same number of elements as alpha has rows (or columns)
    if (size(beta) /= size(alpha, 1)) then
        print *, 'ERROR: alpha does not have enough elements'
        stop
    endif

    ! initialize alpha and beta
    alpha = 0.0_dp
    beta = 0.0_dp

    ! initialize sizes
    i_size = size(alpha,2)
    j_size = size(alpha,1)

    ! initializing the number of data and parameters
    n_data = size(uncertainties)

    ! do loop from notes but in terms of alpha and beta(b-vector)
    do k=1,n_data
        call calculate_linear_terms(n_protons(k), n_neutrons(k), linear_terms)
        do i=1, i_size
            beta(i) = (linear_terms(i) * exp_values(k) / uncertainties(k)**2) + beta(i)
            do j=1, j_size
                alpha(i,j) = (linear_terms(i) * linear_terms(j) / uncertainties(k)**2) + alpha(i,j)
            enddo
        enddo
    enddo
end subroutine construct_alpha_beta

!-----------------------------------------------------------------------
!! Subroutine: calculate_linear_termns
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine is responsible for determining the linear terms  for 
!! the alpha and beta matrices. Each term is equal to a corresponding term
!! given by the semi-empirical mass formula. So the first linear term
!! is equal to the volume_term, whereas the second linear term is equal to
!! the surface term, and so on.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                integer     number of protons in an isotope
!! N                integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! linear_terms        real        Array containing the linear terms in the semi-empirical mass formula
!-----------------------------------------------------------------------
subroutine calculate_linear_terms(Z, N, linear_terms)
    implicit none
    integer, intent(in) :: Z, N
    real(dp), intent(out) :: linear_terms(:)

    ! We could write down all the formulas for each term here. However, in
    ! order to keep the code readable and easy to understand  we'll  separate
    ! them into different functions
    linear_terms(1) = volume_term(Z,N)
    linear_terms(2) = surface_term(Z,N)
    linear_terms(3) = asymmetry_term(Z,N)
    linear_terms(4) = coulomb_term(Z,N)
    linear_terms(5) = pairing_term(Z,N)
    ! If you do the extra credit you would add your new term here
    ! linear_terms(6) = my_extra_term(Z,N)
end subroutine calculate_linear_terms

!-----------------------------------------------------------------------
!! function: volume_term
!-----------------------------------------------------------------------
!! By Benjamin Pieczynski
!!
!! This function determines the volume linear term, where c_vol = A,
!! where A = Z + N.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function volume_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N 

    r = Z + N
end function volume_term

!-----------------------------------------------------------------------
!! function: surface_term
!-----------------------------------------------------------------------
!! By Benjamin Pieczynski
!!
!! This function determines the surface linear term, where c_surf = A^2/3
!! where A = Z + N.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        surface term
!-----------------------------------------------------------------------
real(dp) function surface_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    integer :: A

    A = Z + N
    r = (A**(2.0_dp/3.0_dp))
end function surface_term

!-----------------------------------------------------------------------
!! function: asymmetry_term
!-----------------------------------------------------------------------
!! By Benjamin Pieczynski
!!
!! This function determines the asymmetry linear term, where 
!! c_asym = ((N - Z)^2) / A, where A = Z + N.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        asymmetry term
!-----------------------------------------------------------------------
real(dp) function asymmetry_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    integer :: A

    A = Z + N
    r = ((N - Z)**2.0_dp) / A
end function asymmetry_term

!-----------------------------------------------------------------------
!! function: coulomb_term
!-----------------------------------------------------------------------
!! By Benjamin Pieczynski
!!
!! This function determines the coulomb linear term, where 
!! c_coul = (Z(Z - 1)) / A^1/3, where A = Z + N.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        coulomb term
!-----------------------------------------------------------------------
real(dp) function coulomb_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    integer :: A

    A = Z + N
    r = (Z*(Z-1)) / (A**(1.0_dp/3.0_dp))
end function coulomb_term

!-----------------------------------------------------------------------
!! function: pairing_term
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function determines the pairing linear term, where 
!! c_pair = A^-3/4 * delta_ZN, where A = Z + N.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        pairing term
!-----------------------------------------------------------------------
real(dp) function pairing_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    integer :: A
    integer :: delta_ZN = 0

   A = Z + N
    delta_ZN = initialize_delta_ZN(Z, N)
    r = delta_ZN * (A**(-1.0_dp * 3.0_dp / 4.0_dp))

end function pairing_term

!-----------------------------------------------------------------------
!! function: init_delta_ZN
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function determine's the influence of the pairing energy. This is 
!! determined by the protons and neutrons (Z & N). If they are both even
!! then delta_ZN = 1, if they are both odd then they are -1, and if only
!! one is odd then delta_ZN = 0.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! delta_ZN            real        pairing energy influence
!-----------------------------------------------------------------------
real(dp) function initialize_delta_ZN(Z, N) result(delta_ZN)
    implicit none
    integer, intent(in) :: Z, N

    if ((MOD(Z,2) == 0 ) .and. (MOD(N,2) == 0)) then
        delta_ZN = 1
    else if (MOD(Z+N,2) == 1) then
        delta_ZN = 0
    else
        delta_ZN = -1
    endif



end function initialize_delta_ZN
!-----------------------------------------------------------------------
!! Subroutine: print_best_parameters
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine determines the uncertainty for the best parameters and 
!! then prints the best parameter and its uncertainty out onto one line.
!! This happens for every parameter.
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! covariance       real        Array containing covariance matrix
!-----------------------------------------------------------------------
subroutine print_best_parameters(c_parameters, covariance)
    implicit none
    real(dp), intent(in) :: c_parameters(:), covariance(:,:)
    real(dp), allocatable :: sigma(:)
    integer :: i, j
    

    ! This is an example of how to define and use formats
    allocate(sigma(size(c_parameters)))

    sigma = 0.0_dp
    do j = 1, size(c_parameters)
        do i = 1, size(c_parameters)
            sigma(j) = sigma(j) + c_parameters(i) * covariance(j,i)
        enddo
    enddo
    ! How can you use the error formula in the README to calculate the error
    ! bar in each parameter?
       
    print *
    print *, 'Best fit values:              value                 uncertainty'
    print 1, ' Volume parameter:   ', c_parameters(1), sigma(1)
    print 1, ' Surface parameter:  ', c_parameters(2), sigma(2)
    print 1, ' Asymmetry parameter:', c_parameters(3), sigma(3)
    print 1, ' Coulomb parameter:  ', c_parameters(4), sigma(4)
    print 1, ' Pairing term:       ', c_parameters(5), sigma(5)
        ! add a 6th term for the additional credit portion

1 format(a,f15.8,e28.16)
end subroutine print_best_parameters



!-----------------------------------------------------------------------
!! function: semi_empirical_mass
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function determines the theoretical binding energy for each 
!! proton and neutron pair given our c terms(c_vol, c_surf,... etc).
!! This is based in the semi-empirical mass formula hence each term is added
!! together after the linear term is multiplied by the correct c term.
!!----------------------------------------------------------------------
!! Input:
!!
!! c    real        Array containing the parameters of the semi-empirical mass formula
!! Z    integer     number of protons in an isotope
!! N    integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r    real        Binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_mass(c, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: c(:)
    integer, intent(in) :: Z, N
    real(dp), allocatable :: linear_terms(:)

    ! we need to allocate for linear_terms
    allocate(linear_terms(1:size(c)))
    ! You can call the calculate_linear_termns subroutine and use its output
    ! to calculate the binding energy
    call calculate_linear_terms(Z, N, linear_terms)
    r = c(1) * linear_terms(1) + c(2) * linear_terms(2) &
    + c(3) * linear_terms(3) + c(4) * linear_terms(4) &
     + c(5) * linear_terms(5) !+ c(6) * linear_terms(6)
end function semi_empirical_mass

!-----------------------------------------------------------------------
!! function: semi_empirical_error
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function determines the theoretical binding energy error for each 
!! proton and neutron pair given our c terms(c_vol, c_surf,... etc).
!! This is can be determined by sigma^2 =  sum of the linear_terms(i) * linear_terms(j)
!! * covariance(i,j). We then take the square root to find the error.
!!----------------------------------------------------------------------
!! Input:
!!
!! covariance   real        2D array containing the parameters' covariance matrix
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        statistical uncertainty in the binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_error(covariance, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: covariance(:,:)
    real(dp) :: temp_var
    real(dp), allocatable :: linear_terms(:)
    integer :: i, j, M_size , N_size
    integer, intent(in) :: Z, N


    ! Do you need new functions or subroutine to calculate the derivatives of
    ! the semi-empirical mass formula with respect of each of the c
    ! coefficients or is it already coded? NO

    ! initializing variables
    temp_var = 0.0_dp
    M_size = size(covariance, 2)
    N_size = size(covariance, 1)

    ! allocate for linear terms
    allocate(linear_terms(1:M_size))
    ! call linear terms
    call calculate_linear_terms(Z, N, linear_terms)

    do j=1, M_size
        do i=1, N_size
            temp_var = temp_var + (linear_terms(i) * linear_terms(j) * covariance(i,j))
        enddo
    enddo

    ! resulting error is the square root of the temporary sum/variable
    r = sqrt(temp_var)
    ! Follow the error formula in the README to calculate the theoretical
    ! error in the binding energy. The formula in the README gives the square
    ! of the error, don't forget to take the square root
end function semi_empirical_error

!-----------------------------------------------------------------------
!! Subroutine: vally_of_stability
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine determines the most stable isotope between a Z of 1 to
!! 118. To do this the program must determine the energy per nucleon for
!! each isotope. This subroutine determines this from another function.
!! The subroutine then loops over each stable isotope to find the most stable
!! which occurs when the energy per nucleon(N+1) is > than energy per nucleon(N-1).
!! These values are then submitted to the array containing the most stable
!! isotope, N_STABLE.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! protons          integer     Array containing the protons
!!
!!----------------------------------------------------------------------
!! Output:
!! 
!! N_STABLE         integer     Array containing the most stable isotopes
!-----------------------------------------------------------------------

subroutine valley_of_stability(c_parameters, protons, N_STABLE)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, allocatable, intent(out) :: N_STABLE(:), protons(:)
    real(dp) :: Ep_Nuc, Ep_Nuc_prev
    integer :: i, ipro, N

    ! allocating the Number of stable isotopes
    allocate(N_STABLE(1:118))
    allocate(protons(1:118))

    ! creating array for protons
    protons = 0
    do ipro = 1, 118
        protons(ipro) = ipro
    enddo

    ! initialize N
    N = 1

    ! Searching for N_STABLE (before Ep_Nuc > Ep_Nuc prev ie. last positve energy)
    do i = 1, 118
        Ep_Nuc_prev = energy_per_nucleon(c_parameters, i, N)
        do  
            Ep_Nuc = energy_per_nucleon(c_parameters, i, N+1)
            if (Ep_Nuc > Ep_Nuc_prev) exit
            N = N + 1
            Ep_Nuc_prev = Ep_Nuc
        end do
        N_STABLE(i) = N
    enddo

    

end subroutine valley_of_stability

!-----------------------------------------------------------------------
!! Function: energy_per_nucleon
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This function determines the energy per nucleon for each binding energy
!! for the stable isotopes. This is computed by dividing the binding 
!! energy by A, where A = Z + N.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! Z                integer     Array containing the number of protons
!! N                integer     Array containing the number of neutrons in an isotope
!!
!!----------------------------------------------------------------------
!! Output:
!! 
!! Ep_Nuc         real          Energy per Nucleon     
!-----------------------------------------------------------------------

real(dp) function energy_per_nucleon(c_parameters, Z, N) result(Ep_Nuc)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(in) :: Z, N 
    real(dp) :: BE 
    integer :: A

    ! BE calculation
    BE = semi_empirical_mass(c_parameters, Z, N)

    ! calculating Ep_Nuc
    A = Z + N
    Ep_Nuc = BE / A

end function energy_per_nucleon

!-----------------------------------------------------------------------
!! Subroutine: theoretical_driplines
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This subroutine finds the location of the neutron driplines. To do this
!! the program calculates the Binding energy using the semi-empirical mass function.
!! The program then uses a do loop over all the stable isotopes (1, 118)
!! and compares the BE(N+1) and BE(N) to search for the drip-line position.
!! We determine the dripline position when the difference between BE(N+1)
!! and BE(N) is less than zero.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!!
!!----------------------------------------------------------------------
!! Output:
!! 
!! N_DRIP        real          drip-line positions     
!-----------------------------------------------------------------------

subroutine theoretical_driplines(c_parameters, N_DRIP)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, allocatable, intent(out) :: N_DRIP(:)
    real(dp) :: BE_prev, BE, Sn
    integer :: N, i

    ! allocate
    allocate(N_DRIP(1:118))

    N = 1

    ! 118 maximum amount of protons in a nuleus
    do i = 1, 118
        BE_prev = semi_empirical_mass(c_parameters, i, N)
        do
            BE = semi_empirical_mass(c_parameters, i, N + 1)
            Sn = BE_prev - BE
            if (Sn < 0) exit
            N = N + 1
            BE_prev = BE 
        enddo
        N_DRIP(i) = N
    enddo

end subroutine theoretical_driplines


end module nuclear_model