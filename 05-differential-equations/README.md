# Orbital Mechanics, 3-Body System

**By: Benjamin Pieczynski | Completed: 11/17/2021**

---

**PROGRAM DESCRIPTION**

This program models the orbits of 2 planets around a primary mass (star) for the duration of a specified final time. To model the orbits, this program uses
a 4th order Runge-Kutta method to solve eight coupled first-order ordinary differential equations derived from Newton's 2nd law (equations of motion). The 
program can use input values from a namelist file (see namelist format), or the program will run a set of default values that can be changed inside of 
read_write.f90. The program will write positional, energy, and angular momentum data for the system during each timestep timestep.

**Namelist Example**

``` namelist
&masses
    primary_mass = 1
    planet_mass_1 = 0.1
    planet_mass_2 = 0.05
/

&initial_conditions
    initial_pos_1 = 0.64, 0.92
    initial_pos_2 = 0.62, 0.1
    initial_vel_1 = 0.42, 0.2
    initial_vel_2 = 0.81, 0.2
/

&solution_parameters
    final_time = 61
    n_steps = 1500
/

&output
    output_file = 'nml_result_1.dat'
/
```

**Equations of Motion**

We consider the motion of two “planets”, with masses <img src="/05-differential-equations/tex/0429e3dd940669f4c728ca27fe915301.svg?invert_in_darkmode&sanitize=true" align=middle width=20.985647099999987pt height=14.15524440000002pt/> and <img src="/tex/d9ad343d20544ab9321998ec5d49eba3.svg?invert_in_darkmode&sanitize=true" align=middle width=20.985647099999987pt height=14.15524440000002pt/> and
positions <img src="/05-differential-equations/tex/362188493fdc1549d8b84ac63febac55.svg?invert_in_darkmode&sanitize=true" align=middle width=52.29465614999999pt height=24.65753399999998pt/> and <img src="/05-differential-equations/tex/03156071bd8c54e87abf6d36d80a39c5.svg?invert_in_darkmode&sanitize=true" align=middle width=52.29465614999999pt height=24.65753399999998pt/> about a stationary primary of mass <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/>
located at the origin. All three objects interact through gravitational
attraction, although the primary does not move. For simplicity, we take
Newton’s gravitational constant <img src="/05-differential-equations/tex/4de6f2f08ce4e29f97f32d2fb0410b38.svg?invert_in_darkmode&sanitize=true" align=middle width=43.06148384999999pt height=22.465723500000017pt/>. From Newton’s second law, we have:

<p align="center"><img src="/05-differential-equations/tex/80a859a0c64d64c888e9c5aaed173dbf.svg?invert_in_darkmode&sanitize=true" align=middle width=322.0095945pt height=180.3696015pt/></p>

<!-- <p align="center"><img src="/05-differential-equations/tex/e59b9bade1ace753c7636399d09c07ae.svg?invert_in_darkmode&sanitize=true" align=middle width=319.26996255pt height=40.160877899999996pt/></p>
<p align="center"><img src="/05-differential-equations/tex/54f4ffea2052a54f4b22ecca82eed80e.svg?invert_in_darkmode&sanitize=true" align=middle width=313.92773609999995pt height=40.160877899999996pt/></p>
<p align="center"><img src="/05-differential-equations/tex/bd1783500a3fb0b0d9ddc07ee6d53a24.svg?invert_in_darkmode&sanitize=true" align=middle width=319.26996255pt height=40.160877899999996pt/></p>
<p align="center"><img src="/05-differential-equations/tex/3427051ea2c055efbbacf3981fa4c376.svg?invert_in_darkmode&sanitize=true" align=middle width=313.92773609999995pt height=40.160877899999996pt/></p>
 -->
where 

<p align="center"><img src="/05-differential-equations/tex/9e67799c4be7f5a336882bfe3f60deb5.svg?invert_in_darkmode&sanitize=true" align=middle width=517.12669965pt height=29.58934275pt/></p>

In order to make this suitable for Runge Kutta integration, we introduce
velocities and also divide out masses where we can:

<p align="center"><img src="/05-differential-equations/tex/0e11800ea6e92ee24d0ca2b3a0f86f5a.svg?invert_in_darkmode&sanitize=true" align=middle width=155.33769405pt height=74.19953475pt/></p>

<p align="center"><img src="/05-differential-equations/tex/d58f56bd4c2b0353c9d63d2b2bba8ff2.svg?invert_in_darkmode&sanitize=true" align=middle width=253.26728625pt height=172.50820785pt/></p>


We can also write the total energy, which ought to be conserved:

<p align="center"><img src="/05-differential-equations/tex/4468fa7d1820d0dfe481aae7d554911c.svg?invert_in_darkmode&sanitize=true" align=middle width=510.84627329999995pt height=36.09514755pt/></p>

The total angular momentum of the system is then just the result of the m1*(x1*vy1 - y1*vx1) + m2*(x2*vy2 - y2*vx2).

**4th Order Runge Kutta**

Given the following system of coupled, first-order differential equations:

<p align="center"><img src="/05-differential-equations/tex/ce8bb26ba60e55710b3349a98036a20d.svg?invert_in_darkmode&sanitize=true" align=middle width=108.80473725pt height=33.81208709999999pt/></p>

where <img src="/05-differential-equations/tex/f271c5c5eb7c297ef173986da8101744.svg?invert_in_darkmode&sanitize=true" align=middle width=26.68950404999999pt height=24.65753399999998pt/> is a ***vector*** containing the generalized dependent
variables (which can be positions, velocities or any other type of variables),
and <img src="/05-differential-equations/tex/9ee0b50c5208d83b6078b1dfbffdc738.svg?invert_in_darkmode&sanitize=true" align=middle width=13.02241379999999pt height=32.16441360000002pt/> is a vector function that depends on the generalized variables
and time, 4th order Runge Kutta can be implemented with:

<p align="center"><img src="/05-differential-equations/tex/d4f9efd5b355fadef10c19d1ee7af974.svg?invert_in_darkmode&sanitize=true" align=middle width=285.7352256pt height=189.78578685pt/></p>

where <img src="/05-differential-equations/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> is the time interval in your integration. It is vitally important to
complete each line first before going to the next line. Be sure to pay
attention to factors of <img src="/05-differential-equations/tex/47d54de4e337a06266c0e1d22c9b417b.svg?invert_in_darkmode&sanitize=true" align=middle width=6.552545999999997pt height=27.77565449999998pt/>, etc.

**INSTALLATION**
  
To install this repository git clone https://github.com/phys580-fall2021/05-differential-equations-benjamin-pieczynski.git.
Locate the repository in the terminal and change to the src directory. Enter the make command to compile the program.

To run, enter: ./planetary name_of_namelist
