
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MonolithicFsiSolver.H"
#include "gtest/gtest.h"
using namespace tubeflow;

class MonolithicFsiSolverTest : public ::testing::Test {
    protected:
        virtual void SetUp() {
            // Physical settings
            scalar r0 = 0.2;
            scalar a0 = M_PI * r0 * r0;
            scalar u0 = 0.1;
            scalar p0 = 0;
            scalar dt = 0.1;
            int N = 5;
            scalar L = 1;
            scalar T = 10;
            scalar rho = 1.225;
            scalar E = 490;
            scalar h = 1.0e-3;
            scalar cmk = std::sqrt(E * h / (2 * rho * r0));

            fsi = new MonolithicFsiSolver(a0, u0, p0, dt, cmk, N, L, T, rho);
        }

        virtual void TearDown() {
            delete fsi;
        }

        MonolithicFsiSolver * fsi;
};

TEST_F(MonolithicFsiSolverTest, object)
{
    ASSERT_TRUE(true);
}

TEST_F(MonolithicFsiSolverTest, run)
{
    fsi->run();

    scalar tol = 1.0e-5;

    fsi::vector pn(5), un(5), an(5);
    pn << -1.04294792e-06, 1.11864827e-04, 2.24772603e-04, 4.52098904e-04, 7.13923087e-04;
    un << 0.1, 0.09999385, 0.10015137, 0.1003671, 0.10058284;
    an << 0.1256636, 0.12567518, 0.12568677, 0.1257101, 0.12573697;

    for (int i = 0; i < 5; i++) {
        ASSERT_NEAR(fsi->pn(i), pn(i), tol);
        ASSERT_NEAR(fsi->un(i), un(i), tol);
        ASSERT_NEAR(fsi->an(i), an(i), tol);
    }
}
