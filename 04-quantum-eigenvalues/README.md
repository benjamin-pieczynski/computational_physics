# Solving Schrödinger Equation with Eigenvalues

**By: Benjamin Pieczynski | Completed: 10/27/2021**

---

**PROGRAM DESCRIPTION**

This program calculates the solution for the one-dimensional Schrödinger Equation (equation 1) using eigenvalues and eigenvectors. The program will find *BOTH* a 
numerical and analytical solution for the Schrödinger Equation given the different potential energies for the Infinite Well(particle in a box), Harmonic oscillator
, and the Woods-Saxon potential. The solutions are determined by a set of user input values for length, the number of steps, and the Woods-Saxon radius. For each variation in potential energies the program: 1) constructs kinetic and potential energet diagaonal and off-diagonal terms, 2) solves for the eigenfunctions and eigen values using *lapack*, 3) Normalizes the eigen-wavefunctions, 4) returns the 3 lowest energies and their corresponding wave functions, 5) writes the results to an ouput file. The analytic solutions are only solved int the Infinite Well and Harmonic oscillator solutions. The ultimate goal of this program is to find the eigenenergies and eigenfunctions for the Woods-Saxon potential (equation 2) and write them to a file.

eq. 1

<p align="left"><img src="/04-quantum-eigenvalues/tex/054c6a538901fafcf642b0f7b48aa019.svg?invert_in_darkmode&sanitize=true" align=middle width=260.6055045pt height=40.11819404999999pt/></p>

eq. 2

<p align="left"><img src="/04-quantum-eigenvalues/tex/8355156bca5e02c2b790b4a34942921f.svg?invert_in_darkmode&sanitize=true" align=middle width=206.4436275pt height=37.73900955pt/></p>


**NUMERICAL SOLUTION**

To accomplish the numerical portion, the program turns the wavefunction <img src="/04-quantum-eigenvalues/tex/535b15667b86f1b118010d4c218fecb9.svg?invert_in_darkmode&sanitize=true" align=middle width=12.785434199999989pt height=22.465723500000017pt/> into a column vector by discretizing the function on a lattice. This is discritized in a box from -L to L with a step size of `dx = (2L)/(N-1)` and 'x_vector = -L + dx*(i-1)'. For the Hamiltonian the program creates a tridiagonal matrix with the kinetic energy, where the kinetic energy is a constant `+hbar**2/mass/dx**2`, and off-diagonal being a constant `-0.5*hbar**2/mass/dx**2`. The kinetic energy also needs to be multiplied by the following <img src="/04-quantum-eigenvalues/tex/e72fb7bb3f922789338b5ecbdef00cfc.svg?invert_in_darkmode&sanitize=true" align=middle width=94.5721095pt height=33.45973289999998pt/>. The potential energy is entirely located on the diagonal of another matrix and will vary for different potential energies. For the infinite well the potential energy value is zero. For the Harmonic oscillator the value is <img src="/04-quantum-eigenvalues/tex/ad3b56ed065c54407f0808eadaa74d63.svg?invert_in_darkmode&sanitize=true" align=middle width=95.45020154999999pt height=33.45973289999998pt/>. For the Woods Saxon potential the potential energey is described by equation number 2. After constructing the diagonals the program will then solve for the eigen values(energies) and eigenvectors(wavefunctions) by using the an eigensolver function that uses the *lapack* package 'dstev'. These values are then sent to a subroutine to normalize the wavefunction and then other subroutines that will write them to files.

**ANALYTIC SOLUTION**

To find the analytic solutions, the program applys analytic equations for the Infinite well and Harmonic oscillator potentials and then solves them. To solve these functions we hardcode the parameters found in each equation. Since we only want the 3 lowest solutions, the program will select the 3 lowest energies for comparison with the 3 lowest energies from the numerical section. Equation 3 shows analytic solution for the particle in the box. Equation 4 shows the analytic solution for the analytic solution for the harmonic oscillator.

eq. 3
<p align="left"><img src="/04-quantum-eigenvalues/tex/0ce6560583ba3b1af00c79f559b85de8.svg?invert_in_darkmode&sanitize=true" align=middle width=102.28825859999999pt height=35.77743345pt/></p>

eq. 4
<p align="left"><img src="/04-quantum-eigenvalues/tex/46171ec8ca6d02ef277e160714b2f8e3.svg?invert_in_darkmode&sanitize=true" align=middle width=147.85895519999997pt height=39.452455349999994pt/></p>

**OUTPUT FILES LIST**

1. Infinite Well -> 'infinite_well_wf.dat'
2. Harmonic Oscillator -> 'harmonic_oscillator_wf.dat'
3. Woods-Saxon -> 'woods_saxon_wf.dat'
4. Woods-Saxon R 2 through 10 -> 'woods_saxon_ener.dat'

**INSTALLATION**
  
To install this repository git clone https://github.com/phys580-fall2021/04-quantum-eigenvalues-benjamin-pieczynski.git
