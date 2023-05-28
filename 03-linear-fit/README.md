**PROGRAM: nuclear_energies**

**BY: Benjamin Pieczynski**

**COMPLETED: 10/11/2021**

---

**PROGRAM DESCRIPTION**

This program determines the most stable isotopes and neutron drip-line positions. To accomplish this we use the semi-emperical mass formula(eq. 1) to calculate the 
binding energies for different elements and their isotopes. This program has parallels to investigating stellar nucleosynthesis. The experimental values for
every element and isotope has not been directly measured, so we need to devlope a nuclear binding energy model in order to determine the theoretical values for
the isotope's whose binding energy cannot be directly measured. This program uses the existing experimental values and errors for each element and different isotopes
in order to develope a binding energy model that can then be used to determine the stable isotopes (valley of stability) and the nuclear drip-line positions.

eq 1.
<p align = "center">
<a href="https://www.codecogs.com/eqnedit.php?latex=BE(Z,N)&space;=&space;c_{vol}A&space;&plus;&space;c_{surf}A^{2/3}&space;&plus;&space;c_{asym}\frac{(N-Z)^2}{A}&space;&plus;&space;c_{coul}\frac{Z(Z-1)^2}{A^{1/3}}&space;&plus;c_{pair}A^{-3/4}\delta(Z,N)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?BE(Z,N)&space;=&space;c_{vol}A&space;&plus;&space;c_{surf}A^{2/3}&space;&plus;&space;c_{asym}\frac{(N-Z)^2}{A}&space;&plus;&space;c_{coul}\frac{Z(Z-1)^2}{A^{1/3}}&space;&plus;c_{pair}A^{-3/4}\delta(Z,N)" title="BE(Z,N) = c_{vol}A + c_{surf}A^{2/3} + c_{asym}\frac{(N-Z)^2}{A} + c_{coul}\frac{Z(Z-1)^2}{A^{1/3}} +c_{pair}A^{-3/4}\delta(Z,N)" /></a>
</p>

Where <a href="https://www.codecogs.com/eqnedit.php?latex=A&space;=&space;Z&space;&plus;&space;N" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A&space;=&space;Z&space;&plus;&space;N" title="A = Z + N" /></a>
is the mass number, and <a href="https://www.codecogs.com/eqnedit.php?latex=\delta&space;(Z,N)&space;=&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta&space;(Z,N)&space;=&space;1" title="\delta (Z,N) = 1" /></a> 
if Z and N are both even, = -1 if Z and N are both odd and = 0 if A is odd.

**PROGRAM OVERVIEW**

The first step to developing a model for the binding energies is to determine the best fit for each parameter. The parameters the program solves for is the Volume Term,
Surface Term, Symmetry Term, Coloumb Term, and Pairing Term. Before this can be done we need to input a file that contains the experimental values for all known elements
and their isotopes. The program will prompt the user for a for a file with experimental data, I recommend uploading a file with 2016 experimental results EXPERIMENT_AME2016.dat. 
The program will then read in the neutrons, protons, experimental values, and experimental uncertainty associated with each element and its isotope. After acquiring these values
the program solves a linear system of equations to determine the best-fit parameters for each term (eq. 2). Specifically this involves inverting the alpha matrix in order
to solve for the x vector. To accomplish this the program also has to create an 5x5 alpha matrix (eq. 3) and beta vector with 5 terms.

eq 2.
<p align = "center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;\cdot&space;\vec{x}&space;=&space;\vec{\beta}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;\cdot&space;\vec{x}&space;=&space;\vec{\beta}" title="\alpha \cdot \vec{x} = \vec{\beta}" /></a>
</p>

eq 3.
<p align = "center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\alpha_{ij}&space;=&space;\sum_{k=1}&space;{f}_{i}(Z_{k},&space;N_{k}){f}_{j}(Z_{k},N_{k})&space;/&space;\sigma^2_{k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha_{ij}&space;=&space;\sum_{k=1}&space;{f}_{i}(Z_{k},&space;N_{k}){f}_{j}(Z_{k},N_{k})&space;/&space;\sigma^2_{k}" title="\alpha_{ij} = \sum_{k=1} {f}_{i}(Z_{k}, N_{k}){f}_{j}(Z_{k},N_{k}) / \sigma^2_{k}" /></a>

After solving the linear terms the program can determine both the drip-line locations and the valley of stability positions, which include the most stable isotopes.
To do this we need to plug in values for each isotope from 1 to 118 in order to check each binding energy. To get the valley of stability positions we need to have the
energy per nucleon which requires dividing the binding energy by A.
  
**INSTALLATION**
  
To install this repository git clone https://github.com/phys580-fall2021/03-linear-fit-benjamin-pieczynski.git
 
