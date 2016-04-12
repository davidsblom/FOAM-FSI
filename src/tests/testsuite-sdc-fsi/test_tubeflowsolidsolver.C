
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowLinearizedSolidSolver.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "SDC.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "RBFCoarsening.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "AitkenPostProcessing.H"
#include "AbsoluteConvergenceMeasure.H"

TEST( SDCLinearizedSolidTest, order )
{
    scalar r0 = 0.2;
    scalar h = 1.0e-3;
    scalar L = 1;
    scalar rho_s = 1.225;
    scalar E0 = 490;
    scalar G = 490;
    scalar nu = 0.5;
    scalar p0 = 1e-7;
    scalar T = 1;

    int N = 10;
    int nbComputations = 3;

    std::deque<std::shared_ptr<tubeflow::TubeFlowLinearizedSolidSolver> > solidSolvers;
    std::deque<int> nbTimeStepsList;

    for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
    {
        int nbTimeSteps = std::pow( 2, iComputation );
        scalar dt = T / nbTimeSteps;

        std::shared_ptr<tubeflow::TubeFlowLinearizedSolidSolver> solid( new tubeflow::TubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0, T ) );
        solid->p.fill( p0 );

        solid->run();

        solidSolvers.push_back( solid );
        nbTimeStepsList.push_back( nbTimeSteps );
    }

    fsi::vector ref = solidSolvers.back()->u;

    std::deque<scalar> errors;

    for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
    {
        fsi::vector data;

        data = solidSolvers.at( iComputation )->u;

        scalar error = (ref - data).norm() / data.norm();
        std::cout << "error = " << error << std::endl;
        errors.push_back( error );
    }

    for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
    {
        scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
        std::cout << "order = " << order << std::endl;
        ASSERT_GE( order, 1 );
    }
}
