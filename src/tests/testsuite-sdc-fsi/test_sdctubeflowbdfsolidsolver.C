
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowBDFLinearizedSolidSolver.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "SDC.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "AitkenPostProcessing.H"
#include "AbsoluteConvergenceMeasure.H"

TEST( SDCFsiBDFSolidTest, order )
{
    scalar r0 = 0.2;
    scalar a0 = M_PI * r0 * r0;
    scalar u0 = 0.1;
    scalar p0 = 0;
    scalar L = 1;
    scalar T = 1;
    scalar rho_f = 1.225;
    scalar rho_s = 1.225;
    scalar E0 = 490;
    scalar G = 490;
    scalar h = 1.0e-3;
    scalar nu = 0.5;
    scalar cmk = std::sqrt( E0 * h / (2 * rho_f * r0) );

    int N = 5;
    bool parallel = false;
    int extrapolation = 0;
    scalar tol = 1.0e-3;
    int maxIter = 50;
    scalar initialRelaxation = 1.0e-3;
    int maxUsedIterations = 50;
    int nbReuse = 0;
    scalar singularityLimit = 1.0e-13;
    int reuseInformationStartingFromTimeIndex = 0;
    bool scaling = false;
    scalar beta = 0.01;
    bool updateJacobian = false;
    int nbComputations = 3;

    for ( int timeOrder : {1, 2} )
    {
        std::deque<std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> > fluidSolvers;
        std::deque<std::shared_ptr<tubeflow::SDCTubeFlowBDFLinearizedSolidSolver> > solidSolvers;
        std::deque<int> nbTimeStepsList;

        for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
        {
            int nbTimeSteps = 4 * std::pow( 2, iComputation );
            std::cout << "nbTimeSteps = " << nbTimeSteps << std::endl;
            scalar dt = T / nbTimeSteps;

            std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
            std::shared_ptr<tubeflow::SDCTubeFlowBDFLinearizedSolidSolver> solid( new tubeflow::SDCTubeFlowBDFLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0, T, timeOrder ) );

            shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, 0, 0 ) );
            shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, 1, 0 ) );

            std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
            convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new ResidualRelativeConvergenceMeasure( 0, true, tol ) ) );
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new AbsoluteConvergenceMeasure( 0, true, 1.0e-14 ) ) );

            shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

            shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

            std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( fluid );
            std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( solid );

            assert( sdcFluidSolver );
            assert( sdcSolidSolver );

            std::shared_ptr<fsi::SDCFsiSolver> fsiSolver( new fsi::SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing, extrapolation ) );

            int nbNodes = 3;

            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
            quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::Uniform<scalar>( nbNodes ) );

            std::shared_ptr<sdc::SDC> sdc( new sdc::SDC( fsiSolver, quadrature, 1.0e-12, 10, 50 ) );

            sdc->run();
            ASSERT_TRUE( sdc->isConverged() );

            fluidSolvers.push_back( fluid );
            solidSolvers.push_back( solid );
            nbTimeStepsList.push_back( nbTimeSteps );
        }

        for ( int i = 0; i < 2; i++ )
        {
            fsi::vector ref;

            if ( i == 0 )
                ref = solidSolvers.back()->r;
            else
                ref = solidSolvers.back()->u;

            std::deque<scalar> errors;

            for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
            {
                fsi::vector data;

                if ( i == 0 )
                    data = solidSolvers.at( iComputation )->r;
                else
                    data = solidSolvers.at( iComputation )->u;

                scalar error = (ref - data).norm() / data.norm();
                std::cout << "error = " << error << std::endl;
                errors.push_back( error );
            }

            for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
            {
                scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
                std::cout << "solid order = " << order << std::endl;
                ASSERT_GE( order, timeOrder );
            }
        }

        for ( int i = 0; i < 2; i++ )
        {
            fsi::vector ref;

            if ( i == 0 )
                ref = fluidSolvers.back()->u;
            else
                ref = fluidSolvers.back()->a;

            std::deque<scalar> errors;

            for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
            {
                fsi::vector data;

                if ( i == 0 )
                    data = fluidSolvers.at( iComputation )->u;
                else
                    data = fluidSolvers.at( iComputation )->a;

                scalar error = (ref - data).norm() / data.norm();
                std::cout << "error = " << error << std::endl;
                errors.push_back( error );
            }

            for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
            {
                scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
                std::cout << "fluid order = " << order << std::endl;
                ASSERT_GE( order, 2.99 );
            }
        }
    }
}
