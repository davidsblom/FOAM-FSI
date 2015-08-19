/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    pimpleDyMFoam.C

Description
    Transient solver for incompressible, flow of Newtonian fluids
    on a moving mesh using the PIMPLE (merged PISO-SIMPLE) algorithm.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <string>
#include <memory>
#include <ctime>
#include <math.h>
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // == options == //
    int Nx=1e5;
    int Ny=200;

    if(argc==2)
    {
        Info << "Not enough arguments for mesh. Specify at least 2 (Nx, Ny)" << endl;
    }
    else if(argc>2)
    {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        Info << "Nx, Ny = " << Nx << ", " << Ny << endl;
    }
    int N=Nx*Ny;

    Info << "Building position matrix" << endl;
    Eigen::MatrixXd x1(Ny,3);
    Eigen::MatrixXd M1(Nx,Ny);
    Eigen::MatrixXd M2(int(Nx/2),Ny);
    Eigen::MatrixXd M4(int(Nx/4),Ny);
    Eigen::MatrixXd b1(Nx,3);
    Eigen::MatrixXd b2(int(Nx/2),3);
    Eigen::MatrixXd b4(int(Nx/4),3);

    //Standard testing
    std::clock_t t = std::clock();
    b1 = M1*x1;
    t = std::clock() -t;
    double t1 = 1.0*t/CLOCKS_PER_SEC;
    t = std::clock();
    b2 = M2*x1;
    t = std::clock() -t;
    double t2 = 1.0*t/CLOCKS_PER_SEC;
    t = std::clock();
    b4 = M4*x1;
    t = std::clock() -t;
    double t4 = 1.0*t/CLOCKS_PER_SEC;
    t = std::clock();

    std::cout << "Testing matrix X vector of " << Nx << " by " << Ny << nl;
    std::cout << "t1 | t2 | t4 = " << t1 << " | " << t2 << " | " << t4 << nl << nl;

    //No alias testing
    b1.noalias() = M1*x1;
    t = std::clock() -t;
    double t1na = 1.0*t/CLOCKS_PER_SEC;
    t = std::clock();
    b2.noalias() = M2*x1;
    t = std::clock() -t;
    double t2na = 1.0*t/CLOCKS_PER_SEC;
    t = std::clock();
    b4.noalias() = M4*x1;
    t = std::clock() -t;
    double t4na = 1.0*t/CLOCKS_PER_SEC;
    t = std::clock();

    std::cout << "Testing matrix X vector of " << Nx << " by " << Ny << " by " << 3 <<  " with noalias" << nl;
    std::cout << "t1 | t2 | t4 = " << t1na << " | " << t2na << " | " << t4na << nl << nl;


    return 0;
}


// ************************************************************************* //
