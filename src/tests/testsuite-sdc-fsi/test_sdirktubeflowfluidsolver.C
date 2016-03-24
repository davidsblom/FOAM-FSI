
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "GaussRadau.H"
#include "SDC.H"
#include "ESDIRK.H"
#include "AdaptiveTimeStepper.H"

class SDIRKFluidSolverTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            scalar r0 = 0.2;
            scalar a0 = M_PI * r0 * r0;
            scalar u0 = 0.1;
            scalar p0 = 0.1;
            scalar dt = 0.1;
            int N = 5;
            scalar L = 1;
            scalar T = 1;
            scalar dx = L / N;
            scalar rho = 1.225;
            scalar E = 490;
            scalar h = 1.0e-3;
            scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
            scalar tau = u0 * dt / L;

            ASSERT_NEAR( tau, 0.01, 1.0e-13 );
            ASSERT_TRUE( dx > 0 );

            std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );

            std::shared_ptr<sdc::AdaptiveTimeStepper>adaptiveTimeStepper( new sdc::AdaptiveTimeStepper( false ) );
            esdirk = std::shared_ptr<sdc::ESDIRK> ( new sdc::ESDIRK( fluid, "SDIRK3", adaptiveTimeStepper ) );
        }

        virtual void TearDown()
        {
            esdirk.reset();
        }

        std::shared_ptr<sdc::ESDIRK> esdirk;
};

TEST_F( SDIRKFluidSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( SDIRKFluidSolverTest, timeStep )
{
    esdirk->solveTimeStep( 0 );
}

TEST_F( SDIRKFluidSolverTest, run )
{
    esdirk->run();
}

TEST( SDIRKFluidTest, order )
{
    scalar r0 = 0.2;
    scalar a0 = M_PI * r0 * r0;
    scalar u0 = 0.1;
    scalar p0 = 0;
    int N = 5;
    scalar L = 1;
    scalar T = 1;
    scalar dx = L / N;
    scalar rho = 1.225;
    scalar E = 490;
    scalar h = 1.0e-3;
    scalar cmk = std::sqrt( E * h / (2 * rho * r0) );

    ASSERT_TRUE( dx > 0 );

    int nbComputations = 6;

    std::deque<std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> > fluidSolvers;
    std::deque<int> nbTimeStepsList;

    for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
    {
        int nbTimeSteps = std::pow( 2, iComputation );
        scalar dt = T / nbTimeSteps;

        std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );

        std::shared_ptr<sdc::AdaptiveTimeStepper>adaptiveTimeStepper( new sdc::AdaptiveTimeStepper( false ) );
        std::shared_ptr<sdc::ESDIRK> esdirk( new sdc::ESDIRK( fluid, "ESDIRK53PR", adaptiveTimeStepper ) );

        esdirk->run();

        fluidSolvers.push_back( fluid );
        nbTimeStepsList.push_back( nbTimeSteps );
    }

    fsi::vector pref = fluidSolvers.back()->p;
    std::deque<scalar> errors;

    for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
    {
        scalar error = (pref - fluidSolvers.at( iComputation )->p).norm() / pref.norm();
        errors.push_back( error );
    }

    for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
    {
        scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
        std::cout << "order = " << order << std::endl;
        ASSERT_GE( order, 3 );
    }
}
