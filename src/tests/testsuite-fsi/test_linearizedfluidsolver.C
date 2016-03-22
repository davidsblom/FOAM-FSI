
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearizedFluidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;

class TubeFlowLinearizedFluidSolverTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            scalar p0 = 0;
            scalar r0 = 0.2;
            scalar u0 = 0.1;
            scalar dt = 0.1;
            scalar rho = 1.225;
            scalar L = 1;
            scalar T = 1;
            scalar h = 1.0e-3;
            scalar E = 490;
            int N = 5;

            fluid = new tubeflow::TubeFlowLinearizedFluidSolver( N, p0, r0, u0, rho, E, h, T, dt, L );
        }

        virtual void TearDown()
        {
            delete fluid;
            fluid = NULL;
        }

        TubeFlowLinearizedFluidSolver * fluid;
};

TEST_F( TubeFlowLinearizedFluidSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( TubeFlowLinearizedFluidSolverTest, solve )
{
    fluid->initTimeStep();

    // Initialize variables
    int N = 5;
    fsi::vector p( N ), a( N );
    p.setZero();
    scalar r0 = 0.2;
    scalar a0 = M_PI * r0 * r0;
    a.fill( a0 );

    fluid->solve( a, p );
    fluid->finalizeTimeStep();
}
