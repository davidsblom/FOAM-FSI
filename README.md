FOAM-FSI
========

Features
-----------

* Partitioned Fluid-Structure Interaction solvers build with standard foam-extend fluid and solid solvers.
* Algorithms for strongly coupled FSI: Aitken under relaxation, IQN-ILS or Anderson acceleration, manifold mapping acceleration, and space mapping acceleration.
* Consistent first and second order time integration schemes for incompressible flow, structural dynamics and FSI.
* Efficient implementation of mesh motion solver based on radial basis function interpolation.
* Run fluid-structure-acoustics simulations with the coupling tool preCICE.
* Block-coupled transient incompressible flow solver.

Installation
-----------

First, compile foam-extend-3.1, the nextRelease branch.

To compile the FSI library:

``` bash
git clone https://github.com/davidsblom/FOAM-FSI.git
cd FOAM-FSI
cd src/thirdParty
./Allmake
cd ../../
./Allwmake
```

Prerequisites
-----------

gcc 4.8 or higher due to C++11 features.

Test suite
-----------

A test suite is available implementing a large number of unit and functional tests within the googletest framework. In case you run into problems with the FSI solver, preferably run the test suite and attach a log of the test results to your bug report.
To run the test suite:
``` bash
tests
```

Credits
-----------

The FOAM-FSI project has been started by David Blom. Many other people have also contributed. Here is an alphabetical list:

Philip Cardiff  
Thijs Gillebaart
