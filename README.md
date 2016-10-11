FOAM-FSI
========

[![wercker status](https://app.wercker.com/status/80ce1a6b7a23e394c5f07a3f27b5b43e/m/master "wercker status")](https://app.wercker.com/project/bykey/80ce1a6b7a23e394c5f07a3f27b5b43e)

[![Travis CI status](https://travis-ci.org/davidsblom/FOAM-FSI.svg?branch=master)](https://travis-ci.org/davidsblom/FOAM-FSI)

Features
-----------

*   Partitioned Fluid-Structure Interaction solvers build with standard
    foam-extend fluid and solid solvers.

*   Algorithms for strongly coupled FSI: Aitken under relaxation, IQN-ILS or
    Anderson acceleration, manifold mapping acceleration, and space mapping
    acceleration.

*   Consistent first and second order time integration schemes for
    incompressible flow, structural dynamics and FSI.

*   Efficient implementation of mesh motion solver based on radial basis
    function interpolation.

*   Run fluid-structure-acoustics simulations with the coupling tool preCICE.

*   Block-coupled transient incompressible flow solver.

*   Arbitrary high order time discretization for incompressible flow, solid
    mechanics and FSI simulations.

*   Adaptive radial basis mesh deformation based on coarsening of the
    boundary points.

*   Tested on high performance computing systems, i.e. the SuperMUC Petascale
    System located in Munich. Scalability results show good scaling until
    approximately 2000 MPI ranks.

Installation
-----------

First, compile `foam-extend-3.2`. A detailed description for the compilation of
`foam-extend-3.2` can be found here: <https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-3.2>.
Clone the repository as follows:

``` bash
git clone git://git.code.sf.net/p/foam-extend/foam-extend-3.2
```

It is highly advised to the use the system installed gcc (>= 4.8) and OpenMPI
library for foam-extend and FOAM-FSI. This is due to the fact that the bundled
OpenMPI version of foam-extend disables the fortran compiler, and a set of MPI
routines which are needed when running a fluid-structure-acoustics simulation
with preCICE.

On Ubuntu, the following environment variables should be set before compiling
`foam-extend-3.2`:

``` bash
# Use the system installed libraries
export PARAVIEW_SYSTEM=1
export CUDA_IGNORE=1
export SWAK4FOAM_SYSTEM=1
export WM_MPLIB=SYSTEMOPENMPI
export OPENMPI_DIR=/usr
export OPENMPI_BIN_DIR=$OPENMPI_DIR/bin
```

With this grep command, you can ensure that the `mpicxx` compiler is used
during compilation:

``` bash
sed -i s/"CC          = g++ -m64"/"CC          = mpicxx -m64"/g foam-extend-3.2/wmake/rules/linux64Gcc/c++
```

In case foam-extend is installed in a non-standard location, the environment
variable `$FOAM_INST_DIR` needs to be set to the directory where foam-extend-3.2
is present. So if foam-extend-3.2 is installed at `/home/davidblom/foam-extend-3.2`,
then set `$FOAM_INST_DIR` to `/home/davidblom` like this:

``` bash
export FOAM_INST_DIR=/home/davidblom
```

You can make a distinction between a debug build and a release build. Generally,
it is advised to first setup a test case with a debug build, and then do the
production run with the production build. By default, foam-extend is configured
to build the release version. With the environment variable `WM_COMPILE_OPTION`
you can switch between debug and release builds. So for a debug build:

``` bash
export WM_COMPILE_OPTION=Debug
```

and for a release build:

``` bash
export WM_COMPILE_OPTION=Opt
```

For optimal performance of the production run, it is advised to add the `NDEBUG`
compiler flag in order to disable all assertions as follows:

``` bash
sed -i s/"c++OPT      = -O3"/"c++OPT      = -O3 -DNDEBUG"/g foam-extend-3.2/wmake/rules/linux64Gcc/c++Opt
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

The first time the library is compiled, use the `Allwmake_firstInstall` script
 which installs the necessary settings in `foam-extend-3.2/etc/prefs.sh`, and
 compiles several third party packages. Thereafter, you can use the `Allwmake`
 script.

Prerequisites
-----------

gcc 4.8 or higher due to C++11 features. The gfortran compiler is used for the
compilation of the deal.II library, and is therefore also needed.

MPI: in case the OpenMPI library is used, at least version 1.10 is required
by elemental. Note that foam-extend also needs to be compiled with the same
version.

OS: Currently only Linux based distributions are supported, patches for other
operating systems are welcome.

list of needed packages: scons, openblas (libopenblas-dev on Ubuntu)

Test suite
-----------

Multiple test suite executables are available implementing a large number of unit and functional
tests within the googletest framework. In case you run into problems with the
FSI solver, preferably run the test suite and attach a log of the test results
to your bug report.
The following test suites are available:
* testsuite-dealii
* testsuite-fsi
* testsuite-rbf
* testsuite-sdc
* testsuite-sdc-fsi
* testsuite-spacemapping

To run one set of tests:

``` bash
testsuite-fsi
```

Credits
-----------

The FOAM-FSI project has been started by David Blom. Many other people have also
contributed. Here is an alphabetical list:

Philip Cardiff  
Thijs Gillebaart  
Eline ter Hofstede  
Vahid Kazemi-Kamyab
