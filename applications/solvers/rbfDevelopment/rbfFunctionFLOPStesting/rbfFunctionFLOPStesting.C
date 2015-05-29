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
#include "RBFFunctionInterface.H"
#include "TPSFunction.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Eigen::MatrixXd ematrix;
typedef Eigen::VectorXd evector;

int main(int argc, char *argv[])
{
    // == create "rbfs" == //
    for (int k = 1; k<8; k++)
    {
        Info << nl << " ============================== " << endl;

        int N=1*pow(10,k);
        Info << "N = " << N << endl;

        double R=1.0;
        evector radius(N);
        evector yrbf(N);
        evector dy(N);
        evector dx(N);

        double dN=1.0/N;
        for(int i=1;i<N;i++)
        {
            radius(i)=i*dN;
            dy(i)=std::sin(radius(i));
        }

        std::shared_ptr<rbf::RBFFunctionInterface> rbfFunctionTPS = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        std::shared_ptr<rbf::RBFFunctionInterface> rbfFunctionC2 = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC2Function(R) );

        int t_start = std::clock();
        for ( int i=1; i<N; i++ )
        {
            yrbf(i)=rbfFunctionTPS->evaluate(radius(i));
        }
        int t_end = std::clock();
        double tTPS = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Time to calculate " << N << " TPS function evaluations: " << tTPS << " s" << endl;

        t_start = std::clock();
        for ( int i=1; i<N; i++ )
        {
            yrbf(i)=rbfFunctionC2->evaluate(radius(i));
        }
        t_end = std::clock();
        double tC2 = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Time to calculate " << N << " WendlandC2 function evaluations: " << tC2  << " s"<< endl;

        t_start = std::clock();
        for ( int i=1; i<N; i++ )
        {
            dx(i)=yrbf(i)*dy(i);
        }
        t_end = std::clock();
        double tFLOP = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Time to calculate " << N << " FLOPS: " << tFLOP << " s"<< endl;

        t_start = std::clock();
        dx=yrbf.array()*dy.array();
        t_end = std::clock();
        double tVV = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Time to calculate vector vector product of size " << N << ": " << tVV << " s"<< endl;

        Info << "tTPS/tFLOPS = " << tTPS/tFLOP << endl;
        Info << "tC2/tFLOPS = " << tC2/tFLOP << endl;
        Info << "tTPS/tVV = " << tTPS/tVV << endl;
        Info << "tC2/tVV = " << tC2/tVV << endl;
    }
    Info << " ============================== " << endl;
    return 0;
}


// ************************************************************************* //
