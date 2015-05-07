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
	/*int N=1000;
    rbf::matrix C(N,N);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            C(i,j)=rand();
        }
    }

    int t_start = std::clock();
    //const Eigen::FullPivLU<Eigen::Matrix<double, -1, -1> > Clu = Crbf.fullPivLu();
    rbf::matrix Cinv = C.fullPivLu().inverse();
    int t_end = std::clock();
    double ti = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
    Info << "Inverting C [" << N << "," << N << "]: " << ti << " s" << endl;*/

    // == options == //
    int coarsening=1;
    int Nx=129;
    int Ny=257;
    double R=0.5;
    int nTimes=50;

    // == create "mesh" == //
    //std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC2Function(R) );
    rbf::inverseMatrixLibrary invLib(rbfFunction);
    int N=Nx*Ny;

    double dx=1.0/(Nx-1);
    double dy=1.0/(Ny-1);

    Info << "Building position matrix" << endl;
    rbf::matrix points(N,2);
    rbf::matrix displacement(Nx,2);
    rbf::matrix surfacePoints(Nx,2);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            int k=i+j*Nx;
            points(k,0)=i*dx;
            points(k,1)=j*dy;
        }
        displacement(i,0)=std::rand();
        displacement(i,1)=std::rand();
        surfacePoints(i,0)=i*dx;
        surfacePoints(i,1)=0;
    }

    int Nc=int(Nx/coarsening);
    rbf::matrix controlPoints(Nc,2);
    for(int i=0;i<Nc;i++){
        controlPoints(i,0)=surfacePoints(i*coarsening,0);
        controlPoints(i,0)=surfacePoints(i*coarsening,1);
    }

    Info << "==== mesh properties === " << endl;
    Info << "Ni = " << N << endl;
    Info << "Ns = " << Nx << endl;
    Info << "Nc = " << Nc << endl;

    Info<< nl << "==== Timings - Building ==== " << endl;
    clock_t t_start = std::clock();
    rbf::matrix Crbf(Nc,Nc);
    invLib.getControlMatrix(controlPoints,Crbf);
    clock_t t_end = std::clock();
    double tcb = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
    Info << "Building C [" << Nc << "," << Nc << "]: " << tcb << " s" << endl;

    t_start = std::clock();
    //const Eigen::FullPivLU<Eigen::Matrix<double, -1, -1> > Clu = Crbf.fullPivLu();
    rbf::matrix Cinv = Crbf.fullPivLu().inverse();
    t_end = std::clock();
    double ti = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
    Info << "Inverting C [" << Nc << "," << Nc << "]: " << ti << " s" << endl;

    // ========= Testing matrix calculations ==========//
    /*double totalMatrixVector = 0;
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
    Info << "Symmetric Matrix (lower) vecotr product from 1x1 to " << N << "x" << N << " = " << totalMatrixVectorSymU << " ms" << endl;*/

    return 0;
}


// ************************************************************************* //
