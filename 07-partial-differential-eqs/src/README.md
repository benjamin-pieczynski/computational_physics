# Time dependent Schrödinger Equations in a Potential

**By: Benjamin Pieczynski | Completed: 12/16/2021**

---

**PROGRAM DESCRIPTION**

This program solves the time-dependent Schrödinger equation (equation 1) for a box of size L using the Crank-Nicolson numerical method. By 
introducing a Hamiltonian operator (equation 2), we can rewrite the Schrödinger equation and discretize 
<img src="/07-partial-differential-eqs/svgs/17fef61e769874dc5430b5eefe47024c.svg?invert_in_darkmode&sanitize=true" align=middle width=48.207885pt height=24.6576pt/> as a vector: <img src="/07-partial-differential-eqs/svgs/e404229d6a7a1d30d4f0c2f027d8dc9d.svg?invert_in_darkmode&sanitize=true" align=middle width=112.42407pt height=32.16444pt/>, where <img src="/07-partial-differential-eqs/svgs/68f3ca6428f40a62356634665d6afac6.svg?invert_in_darkmode&sanitize=true" align=middle width=115.697505pt height=34.33782pt/>.
We then discretize the Hamiltonian as a matrix, where we use the Crank-Nicolson method to split up the Hamiltonian to be half implicit and
the other half explicit. The implicit section of the Hamiltonian would be described as imaginary, but we recast it as a purely real problem.
To accomplish this we combine the imaginary and real portions into a large matrix that is of 2N x 2N size, where N is the size of the Wave function
vector. This matrix is refered to in the program as the time evolution matrix (equation 4). This is constructed after reading the input parameters 
that the program will use to sample lattice points (`delta_x = 2*length/(n_points-1)`). This allows us to set the initial wave function to a
Gaussian (equation 5). To combine the two matrices for the time evolution matrix the program utilizes the `linear_algebra` module to invert the
first matrix. The multiplication is then completed using matmul. Then you can evolve you wave function and store the different snapshots in the
`time_wave_function` array. Then the program calculates the following expectation
values for every time step.

1. The normalization <img src="/07-partial-differential-eqs/svgs/babda9c4be225f4d6cb20d799d988a65.svg?invert_in_darkmode&sanitize=true" align=middle width=60.39396pt height=24.65793pt/>, which should be constant.
2. The position <img src="/07-partial-differential-eqs/svgs/049a396d8b1486a6eda9f828adc329bc.svg?invert_in_darkmode&sanitize=true" align=middle width=189.259455pt height=24.65793pt/>
3. The width <img src="/07-partial-differential-eqs/svgs/6641bf1a02c68083fb3b2b8053a93caa.svg?invert_in_darkmode&sanitize=true" align=middle width=127.539885pt height=29.42511pt/>, where <img src="/07-partial-differential-eqs/svgs/899d2babd87c666fec7168fdcdf6c5f4.svg?invert_in_darkmode&sanitize=true" align=middle width=197.082105pt height=26.76201pt/>

Finally the program writes the results into two different files. One file, specified via namelist or default parameters, should contain
the expectation values (one per column) as a function of time. The second file contains a few snapshots of the probability density 
<img src="/07-partial-differential-eqs/svgs/7872e8aeeb0376e1d2b14b69fbc75010.svg?invert_in_darkmode&sanitize=true" align=middle width=30.67944pt height=24.6576pt/>.

For a more advanced version of the code the wavefunction can be put inside a harmonic well , <img src="/07-partial-differential-eqs/svgs/3099cd286b75db93545583aff185dbf8.svg?invert_in_darkmode&sanitize=true" align=middle width=92.860845pt height=27.77577pt/>.
In the basic case the program will initialize the potential to be zero.

**EQUATIONS**

Equation 1)

<p align="left"><img src="/07-partial-differential-eqs/svgs/31c270932cfed8407da7128e5fd05441.svg?invert_in_darkmode&sanitize=true" align=middle width=297.7425pt height=40.118265pt/></p>

Equation 2)

<p align="left"><img src="/07-partial-differential-eqs/svgs/66730f323899e95bbde40b1c5ac99d0d.svg?invert_in_darkmode&sanitize=true" align=middle width=152.922165pt height=33.45969pt/>
  
Equation 3)
 
<p align="left"><img src="/07-partial-differential-eqs/svgs/9e8dd8841293d125d1a409231a2e3bb9.svg?invert_in_darkmode&sanitize=true" align=middle width=181.0611pt height=33.81213pt/></p>

Equation 4)

<p align="left"><img src="/07-partial-differential-eqs/svgs/291beaed1cdda99878eb008def0f7065.svg?invert_in_darkmode&sanitize=true" align=middle width=261.05475pt height=46.07955pt/></p>

Equation 5)

<p align="left"><img src="/07-partial-differential-eqs/svgs/39c380792cdea241cbe0867622a58e39.svg?invert_in_darkmode&sanitize=true" align=middle width=324.76785pt height=40.118265pt/></p>

Equation 6)

<p align="left"><img src="/07-partial-differential-eqs/svgs/291beaed1cdda99878eb008def0f7065.svg?invert_in_darkmode&sanitize=true" align=middle width=261.05475pt height=46.07955pt/></p>


**Namelist Example**

``` namelist
&integration
    length = 10
    n_points = 100
    n_steps = 100
    delta_t = 0.05
/

&wave_function
    width = 0.4
    center = 0.0
/

&oscillator
    k_oscillator = 0.0
/

&output
    time_file = 'time_results.dat'
    density_file = 'density_results.dat'
/
```

**INSTALLATION AND RUNNING THE PROGRAM**
  
To install this repository git clone https://github.com/phys580-fall2021/07-partial-differential-eqs-benjamin-pieczynski.git.
Locate the repository in the terminal and change to the src directory. Enter the make command to compile the program.

To run on a macOS, enter: ./schrodinger name_of_namelist into the terminal

A namelist is included within the src file titled example.namelist
