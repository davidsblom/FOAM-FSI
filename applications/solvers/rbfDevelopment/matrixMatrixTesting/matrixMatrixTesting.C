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
    // == options == //
    int coarsening=4;
    int Nx=256;
    int Ny=256;
    double R=2;
    int nTimes=50;

    if(argc==2)
    {
        Info << "Not enough arguments for mesh. Specify at least 2 (Nx, Ny)" << endl;
    }
    else if(argc==3)
    {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        Info << "Nx, Ny = " << Nx << ", " << Ny << endl;
    }
    else if(argc==4)
    {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        coarsening = atoi(argv[3]);
        Info << "Nx, Ny, coarsening = " << Nx << ", " << Ny << ", " << nTimes << endl;
    }
    else if(argc>4)
    {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        coarsening = atoi(argv[3]);
        nTimes = atoi(argv[4]);
        Info << "Nx, Ny, coarsening, nTimes = " << Nx << ", " << Ny << ", " << coarsening << ", " << nTimes << endl;
    }
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
    rbf::matrix controlDisplacement(Nc,2);
    for(int i=0;i<Nc;i++){
        controlPoints(i,0)=surfacePoints(i*coarsening,0);
        controlPoints(i,1)=surfacePoints(i*coarsening,1);
        controlDisplacement(i,0)=displacement(i,0);
        controlDisplacement(i,1)=displacement(i,1);
    }

    Info << "==== mesh properties === " << endl;
    Info << "Ni = " << N << endl;
    Info << "Ns = " << Nx << endl;
    Info << "Nc = " << Nc << endl;

    // == create "mesh" == //
    for ( int k = 0; k<5;k++)
    {
        std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
        if(k==0)
        {
            Info << nl << "Using TPS rbf function" << endl;
            rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        }
        else if(k==1)
        {
            Info << nl << "Using WendlandC0 rbf function" << endl;
            rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC0Function(R) );
        }
        else if(k==2)
        {
            Info << nl << "Using WendlandC2 rbf function" << endl;
            rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC2Function(R) );
        }
        else if(k==3)
        {
            Info << nl << "Using WendlandC4 rbf function" << endl;
            rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC4Function(R) );
        }
        else
        {
            Info << nl << "Using WendlandC6 rbf function" << endl;
            rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC6Function(R) );
        }

        rbf::inverseMatrixLibrary invLib(rbfFunction);

        Info<< nl << "==== Timings - Building ==== " << endl;
        clock_t t_start = std::clock();
        rbf::matrix Crbf(Nc,Nc);
        invLib.getControlMatrix(controlPoints,Crbf);
        clock_t t_end = std::clock();
        double tcb = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Building C [" << Nc << "," << Nc << "]: " << tcb << " s" << endl;

        t_start = std::clock();
        rbf::matrix Cinv = Crbf.fullPivLu().inverse();
        t_end = std::clock();
        double ti = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Inverting C [" << Nc << "," << Nc << "] using FullPivLU: " << ti << " s" << endl;

        t_start = std::clock();
        rbf::matrix Cinv2 = Crbf.partialPivLu().inverse();
        t_end = std::clock();
        double ti2 = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Inverting C [" << Nc << "," << Nc << "] using PartialPivLU: " << ti2 << " s" << endl;

        t_start = std::clock();
        rbf::matrix H(N,Nc);
        invLib.getEvaluationMatrix(points,controlPoints,H);
        t_end = std::clock();
        double tbb = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Building B [" << N << "," << Nc << "]: " << tbb << " s" << endl;

        t_start = std::clock();
        H*=Cinv;
        t_end = std::clock();
        double tch = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Calculating H [" << N << "," << Nc << "]: " << tch << " s" << endl;
        Info << "Size of H [" << N << "," << Nc << "]: " << N*Nc*8/1000/1000 << " MB" << endl;

        Info<< nl << "==== Timings - Solving ==== " << endl;
        t_start = std::clock();
        for(int i=0;i<nTimes;i++)
        {
            rbf::vector dx=H*controlDisplacement.col(0);
            rbf::vector dy=H*controlDisplacement.col(1);
        }
        t_end = std::clock();
        double tMEM = 1.0*(t_end-t_start)/CLOCKS_PER_SEC;
        Info << "Calculating dx,dy [" << N << "," << "2] "<< nTimes << " times using memory: " << tMEM << " s" << endl;

        double tas=0;double trbfs=0;double tes=0;
        for(int k=0;k<nTimes;k++)
        {
            t_start = std::clock();
            rbf::vector alphax=Cinv*controlDisplacement.col(0);
            rbf::vector alphay=Cinv*controlDisplacement.col(1);
            t_end = std::clock();
            double ta=1.0*(t_end-t_start)/CLOCKS_PER_SEC;

            rbf::vector dx(N);
            rbf::vector dy(N);
            for(int i=0;i<N;i++)
            {
                t_start = std::clock();
                rbf::vector rbfvalue(Nc);
                for(int j=0;j<Nc;j++)
                {
                    double r = std::sqrt(std::pow(points(i,0)-controlPoints(j,0),2) + std::pow(points(i,1)-controlPoints(j,1),2));
                    //double r = ( points.row( i ) - controlPoints.row( j ) ).norm();
                    rbfvalue(j) = rbfFunction->evaluate(r);
                }
                t_end = std::clock();
                double trbf=1.0*(t_end-t_start)/CLOCKS_PER_SEC;

                t_start = std::clock();
                dx(i)=rbfvalue.transpose()*alphax;
                dy(i)=rbfvalue.transpose()*alphay;
                t_end = std::clock();
                double te=1.0*(t_end-t_start)/CLOCKS_PER_SEC;

                trbfs+=trbf;
                tes+=te;
            }
            tas+=ta;
        }
        double tCPU=trbfs+tes;
        Info << "Calculating dx,dy [" << N << "," << "2] "<< nTimes << " times using CPU: " << tCPU << "[" << tas << " , " << trbfs << " , " << tes << "]" << endl;

        /*tas=0;trbfs=0;tes=0;
        for(int k=0;k<nTimes;k++)
        {
            t_start = std::clock();
            rbf::vector alphax=Cinv*controlDisplacement.col(0);
            rbf::vector alphay=Cinv*controlDisplacement.col(1);
            t_end = std::clock();
            double ta=1.0*(t_end-t_start)/CLOCKS_PER_SEC;

            t_start = std::clock();
            rbf::matrix He(N,Nc);
            invLib.getEvaluationMatrix(points,controlPoints,He);
            t_end = std::clock();
            double trbf=1.0*(t_end-t_start)/CLOCKS_PER_SEC;

            t_start = std::clock();
            rbf::vector dx(N);
            rbf::vector dy(N);
            dx=He*alphax;
            dy=He*alphay;
            t_end = std::clock();
            double te=1.0*(t_end-t_start)/CLOCKS_PER_SEC;

            trbfs+=trbf;
            tes+=te;
            tas+=ta;
        }
        double tCPU2=trbfs+tes;
        Info << "Calculating dx,dy [" << N << "," << "2] "<< nTimes << " times using CPU: " << tCPU2 << "[" << tas << " , " << trbfs << " , " << tes << "]" << endl;*/

        Info << "Ratio tCPU/tMEM = " << tCPU/tMEM << endl;
        //Info << "Ratio tCPU2/tMEM = " << tCPU2/tMEM << endl;

    }
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
