# Visualizing Elliptic Curves (and other algebraic groups) over Finite Fields

(Under construction)

This repository contains code designed to help understand elliptic curves over finitel fields
by doing computations involving lattices.

The python code can be found in the ecc/ folder. To use the files, you need to clone the repo and build the wheel.
The imports will trigger an error if you don't build the wheel first, but this can be fixed by deleting the ecc. prefix
from all the module imports. 

The main function is the 'IsogenyClass' class, in the ecc.ecv module. This class is called by specifying a pair of integers
$a,p$, where $p$ is a prime and $a$ is a trace of Frobenius. This will compute an equivalence of categories between the Fp isogeny class, and a class of lattices, encoded as a dictionary. Once this equivalence is computed, one can do several things:
* Compute Mordell-Weil groups as subgroups of an analytic curve, and visualize them using matplotlib
* Compute isogeny graphs of various degrees.

