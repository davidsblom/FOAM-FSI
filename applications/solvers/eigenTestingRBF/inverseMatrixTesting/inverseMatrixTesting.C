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
#include "inverseMatrixLibrary.H"
#include "TPSFunction.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC2Function(1.0) );
    rbf::inverseMatrixLibrary invLib(rbfFunction);

    //unit square with Nx and Ny points
    int Nx=30;
    int Ny=30;
    int N=Nx*Ny;
    Info << "N = " << N << nl << endl;

    double dx=1.0/(Nx-1.0);
    double dy=1.0/(Ny-1.0);

    Info << "Building position matrix" << endl;
    rbf::matrix positions(N,3);
    rbf::matrix displacement(N,1);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            int k=i+j*Nx;
            positions(k,0)=i*dx;
            positions(k,1)=j*dy;
            positions(k,2)=0;
            displacement(k,0)=0.25*(1-std::cos(2*M_PI*i*dx))*(1-std::cos(2*M_PI*j*dy));
        }
    }

    // ========= Testing matrix calculations ==========//
    double totalMatrixVector = 0;
    rbf::matrix Crbf(N,N);
    Info << "Computing rbf control matrix" << endl;
    invLib.getControlMatrix(positions,Crbf);
    Info << "Done computing rbf control matrix" << nl << endl;

    Info << "Starting matrix vector product" << endl;
    for(int ip=1;ip<N;ip++){
        //Start clocking

        rbf::matrix curCrbf(ip,ip);
        rbf::vector curDisp(ip);
        for(int i=0;i<ip;i++){
            for(int j=0;j<ip;j++){
                curCrbf(i,j) = Crbf(i,j);
            }
            curDisp(i) = displacement(i,0);
        }

        rbf::vector result(ip);

        clock_t t_start = std::clock();
        result = curCrbf*curDisp;
        clock_t t_end = std::clock();

        double matrixVector = 1000.0 * (t_end-t_start)/CLOCKS_PER_SEC;
        totalMatrixVector += matrixVector;
    }
    Info << "Matrix vecotr product from 1x1 to " << N << "x" << N << " = " << totalMatrixVector << " ms" << endl;

    // ========= Testing symmetric matrix calculations lower ==========//
    double totalMatrixVectorSymL = 0;
    rbf::matrix CrbfSymL(N,N);
    Info << "Computing rbf control matrix" << endl;
    invLib.getControlMatrixLower(positions,CrbfSymL);
    Info << "Done computing rbf control matrix" << nl << endl;

    Info << "Starting matrix vector product using symmetry" << endl;
    for(int ip=1;ip<N;ip++){
        //Start clocking


        rbf::matrix curCrbf(ip,ip);
        rbf::vector curDisp(ip);
        for(int i=0;i<ip;i++){
            for(int j=0;j<i+1;j++){
                curCrbf(i,j) = CrbfSymL(i,j);
            }
            curDisp(i) = displacement(i,0);
        }


        rbf::vector result(ip);

        clock_t t_start = std::clock();
        //result = curCrbf.selfadjointView<Eigen::Lower>()*curDisp;
        result = curCrbf.triangularView<Eigen::Lower>()*curDisp;
        clock_t t_end = std::clock();

        double matrixVector = 1000.0 * (t_end-t_start)/CLOCKS_PER_SEC;
        totalMatrixVectorSymL += matrixVector;
    }
    Info << "Symmetric Matrix (lower) vecotr product from 1x1 to " << N << "x" << N << " = " << totalMatrixVectorSymL << " ms" << endl;

    // ========= Testing symmetric matrix calculations upper ==========//
    double totalMatrixVectorSymU = 0;
    rbf::matrix CrbfSymU(N,N);
    Info << "Computing rbf control matrix" << endl;
    invLib.getControlMatrixUpper(positions,CrbfSymU);
    Info << "Done computing rbf control matrix" << nl << endl;

    Info << "Starting matrix vector product using symmetry" << endl;
    for(int ip=1;ip<N;ip++){
        //Start clocking


        rbf::matrix curCrbf(ip,ip);
        rbf::vector curDisp(ip);
        for(int i=0;i<ip;i++){
            for(int j=i;j<ip;j++){
                curCrbf(i,j) = CrbfSymU(i,j);
            }
            curDisp(i) = displacement(i,0);
        }


        rbf::vector result(ip);

        clock_t t_start = std::clock();
        //result = curCrbf.selfadjointView<Eigen::Upper>()*curDisp;
        result = curCrbf.triangularView<Eigen::Upper>()*curDisp;
        clock_t t_end = std::clock();

        double matrixVector = 1000.0 * (t_end-t_start)/CLOCKS_PER_SEC;
        totalMatrixVectorSymU += matrixVector;
    }
    Info << "Symmetric Matrix (lower) vecotr product from 1x1 to " << N << "x" << N << " = " << totalMatrixVectorSymU << " ms" << endl;

    /*//Row update
    double totalUpdateInverseMS = 0;
    rbf::matrix oldCinverse(1,1);
    rbf::matrix newpoint(1,3);
    oldCinverse(0,0) = 1.0;//inverse of diagonal value. For compact supported this is 1
    for(int ip=1;ip<N;ip++){
        //Start clocking
        clock_t t_start = std::clock();

        //Make new point
        newpoint(0,0) = positions(ip,0);
        newpoint(0,1) = positions(ip,1);
        newpoint(0,2) = positions(ip,2);

        //Create old positions
        rbf::matrix oldPositions(ip,3);
        for(int j=0;j<ip;j++){
            oldPositions(j,0) = positions(j,0);
            oldPositions(j,1) = positions(j,1);
            oldPositions(j,2) = positions(j,2);
        }

        rbf::matrix inverse(ip+1,ip+1);
        invLib.updateInverse(oldCinverse,newpoint,oldPositions,inverse);

        oldCinverse = inverse;

        clock_t t_end = std::clock();
        double updateInverseMS = 1000.0 * (t_end-t_start)/CLOCKS_PER_SEC;
        totalUpdateInverseMS += updateInverseMS;
    }
    Info << "Updating Inverse from 2x2 to " << N << "x" << N << " = " << totalUpdateInverseMS << " ms" << endl;

    //Full inverse eigen
    double invTimes [N-1];
    double totalInverseEigenMS = 0;
    Info << "Starting eigen inverse (via LU) with eigen" << endl;
    for(int ip=1;ip<N;ip++){
        //Start clocking
        clock_t t_start = std::clock();

        //Preparing function parameters
        rbf::matrix curpositions(ip+1,3);
        for(int j=0;j<ip+1;j++){
            curpositions(j,0)=positions(ip-1,0);
            curpositions(j,1)=positions(ip-1,1);
            curpositions(j,2)=positions(ip-1,2);
        }

        rbf::matrix inverse(ip+1,ip+1);
        //Computing inverse
        invLib.computeInverseEigen(curpositions,inverse);

        //Clocking
        clock_t t_end = std::clock();
        double inverseEigenMS = 1000.0 * (t_end-t_start)/CLOCKS_PER_SEC;
        invTimes[ip-1]=inverseEigenMS;
        totalInverseEigenMS += inverseEigenMS;
    }
    Info << "Time for Eigen inverse of 2x2 to " << N << "x" << N << " = " << totalInverseEigenMS << " ms" << endl;

    //Solve system eigen
    Info << "Starting solve system (via LU) with eigen" << endl;
    double totalSolveEigenMS = 0;
    for(int ip=1;ip<N;ip++){
        //Start clocking
        clock_t t_start = std::clock();

        //Preparing function parameters
        rbf::matrix curpositions(ip+1,3);
        rbf::matrix curdisplacement(ip+1,1);
        for(int j=0;j<ip+1;j++){
            curpositions(j,0)=positions(ip-1,0);
            curpositions(j,1)=positions(ip-1,1);
            curpositions(j,2)=positions(ip-1,2);
            curdisplacement(j,0)=displacement(ip-1,0);
        }
        rbf::matrix solution (ip+1,1);
        invLib.solveSystemEigen(curpositions,curdisplacement,solution);

        clock_t t_end = std::clock();
        double solveEigenMS = 1000.0 * (t_end-t_start)/CLOCKS_PER_SEC;
        totalSolveEigenMS += solveEigenMS;
    }
    Info << "Solving system from 2x2 to " << N << "x" << N << " = " << totalSolveEigenMS << " ms" << endl;*/

    return 0;
}


// ************************************************************************* //
