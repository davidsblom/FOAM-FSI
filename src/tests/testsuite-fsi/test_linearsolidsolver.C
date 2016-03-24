
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;

class TubeFlowLinearSolidSolverTest : public ::testing::Test
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

            scalar dt = 0.1;
            int N = 5;

            solid = new tubeflow::TubeFlowLinearSolidSolver( N, nu, rho, h, L, dt, G, E0, r0 );
        }

        virtual void TearDown()
        {
            delete solid;
            solid = NULL;
        }

        TubeFlowLinearSolidSolver * solid;
};

TEST_F( TubeFlowLinearSolidSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( TubeFlowLinearSolidSolverTest, po )
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
