FOAM-FSI
========

Installation
-----------

First, compile foam-extend-3.1, the nextRelease branch.

To compile the FSI library:

1. cd src/thirdParty
2. ./Allmake
3. cd -
4. ./Allwmake

Prerequisites
-----------

gcc 4.8 or higher due to C++11 features.

Features
-----------

* Partitioned Fluid-Structure Interaction solvers build with standard foam-extend fluid and solid solvers.
* Algorithms for strongly coupled FSI: Aitken under relaxation, IQN-ILS or Anderson acceleration, manifold mapping acceleration, and space mapping acceleration.
* Consistent first and second order time integration schemes for incompressible flow, structural dynamics and FSI.
* Efficient implementation of mesh motion solver based on radial basis function interpolation.

Credits
-----------

The FOAM-FSI project has been started by David Blom. Many other people have also contributed. Here is an alphabetical list:

Philip Cardiff  
Thijs Gillebaart
