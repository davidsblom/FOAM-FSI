FOAM-FSI
========

[![wercker status](https://app.wercker.com/status/80ce1a6b7a23e394c5f07a3f27b5b43e/s/master "wercker status")](https://app.wercker.com/project/bykey/80ce1a6b7a23e394c5f07a3f27b5b43e)
[![codecov.io](https://codecov.io/github/davidsblom/FOAM-FSI/coverage.svg?branch=master)](https://codecov.io/github/davidsblom/FOAM-FSI?branch=master)

Features
-----------

* Partitioned Fluid-Structure Interaction solvers build with standard foam-extend fluid and solid solvers.
* Algorithms for strongly coupled FSI: Aitken under relaxation, IQN-ILS or Anderson acceleration, manifold mapping acceleration, and space mapping acceleration.
* Consistent first and second order time integration schemes for incompressible flow, structural dynamics and FSI.
* Efficient implementation of mesh motion solver based on radial basis function interpolation.
* Run fluid-structure-acoustics simulations with the coupling tool preCICE.
* Block-coupled transient incompressible flow solver.
* Arbitrary high order time discretization for incompressible flow, solid mechanics and FSI simulations.

Installation
-----------

First, compile `foam-extend-3.2`. A detailed description for the compilation of `foam-extend-3.2` can be found here: https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-3.1 which also applies to `foam-extend-3.2`. Note that `foam-extend-3.2` is using the git repository of `foam-extend-3.1`, so use the `v3.2` tag to compile `foam-extend-3.2`:

``` bash
git clone --branch v3.2 git://git.code.sf.net/p/openfoam-extend/foam-extend-3.1 foam-extend-3.2
```

It is highly advised to the use the system installed gcc (>= 4.8) and OpenMPI library for foam-extend and FOAM-FSI. This is due to the fact that the bundled OpenMPI version of foam-extend disables the fortran compiler, and a set of MPI routines which are needed when running a fluid-structure-acoustics simulation with preCICE.

On Ubuntu, the following environment variables should be set before compiling `foam-extend-3.2`:

``` bash
# Use the system installed libraries
export PARAVIEW_SYSTEM=1
export CUDA_IGNORE=1
export SWAK4FOAM_SYSTEM=1
export WM_MPLIB=SYSTEMOPENMPI
export OPENMPI_DIR=/usr
export OPENMPI_BIN_DIR=$OPENMPI_DIR/bin
```

With this grep command, you can ensure that the `mpicxx` compiler is used during compilation:

``` bash
sed -i s/"CC          = g++ -m64"/"CC          = mpicxx -m64"/g foam-extend-3.2/wmake/rules/linux64Gcc/c++
```

Set the environment variable `$FOAM_INST_DIR` in case foam-extend is installed in a non-standard location:

``` bash
export FOAM_INST_DIR=`pwd`
```

Finally, the compilation process of `foam-extend` is started with

``` bash
cd foam-extend-3.2
source etc/bashrc
./Allwmake.firstInstall
```

To compile the `FOAM-FSI` library:

``` bash
git clone --recursive https://github.com/davidsblom/FOAM-FSI.git
cd FOAM-FSI
./Allwmake.firstInstall
```

The first time the library is compiled, use the `Allwmake_firstInstall` script which installs the necessary settings in `foam-extend-3.2/etc/prefs.sh`, and compiles several third party packages. Thereafter, you can use the `Allwmake` script.

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
Eline ter Hofstede  
Vahid Kazemi-Kamyab
