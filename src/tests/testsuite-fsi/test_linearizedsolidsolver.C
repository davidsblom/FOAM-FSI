
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearizedSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;

class TubeFlowLinearizedSolidSolverTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            scalar r0 = 3.0e-3;
            scalar h = 3.0e-4;
            scalar L = 0.126;
            scalar rho = 1000;
            scalar E0 = 4.0e5;
            scalar G = 4.0e5;
            scalar nu = 0.5;
            scalar T = 1;

            scalar dt = 0.1;
            int N = 5;

            solid = new tubeflow::TubeFlowLinearizedSolidSolver( N, nu, rho, h, L, dt, G, E0, r0, T );
        }

        virtual void TearDown()
        {
            delete solid;
            solid = NULL;
        }

        TubeFlowLinearizedSolidSolver * solid;
};

TEST_F( TubeFlowLinearizedSolidSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( TubeFlowLinearizedSolidSolverTest, po )
{
    solid->initTimeStep();

    // Initialize variables
    int N = 5;
    fsi::vector p( N ), a( N );
    p.setZero();
    a.setZero();

    solid->solve( p, a );
    solid->finalizeTimeStep();
}
