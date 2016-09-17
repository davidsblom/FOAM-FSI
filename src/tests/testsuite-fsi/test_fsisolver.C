
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ConvergenceMeasure.H"
#include "FsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using namespace fsi;

class FsiSolverTest : public ::testing::TestWithParam<bool>{
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
            scalar T = 1;
            scalar dx = L / N;
            scalar rho = 1.225;
            scalar E = 490;
            scalar h = 1.0e-3;
            scalar cmk = std::sqrt(E * h / (2 * rho * r0));
            scalar c0 = std::sqrt(cmk * cmk - p0 / (2 * rho));
            scalar kappa = c0 / u0;
            scalar tau = u0 * dt / L;

            // Computational settings
            scalar tol = 1.0e-5;
            int extrapolationOrder = 0;

            bool parallel = GetParam();

            ASSERT_NEAR(tau, 0.01, 1.0e-13);
            ASSERT_NEAR(kappa, 10, 1.0e-13);
            ASSERT_TRUE(dx > 0);

            // Initialize solvers
            shared_ptr<TubeFlowFluidSolver> fluid(new TubeFlowFluidSolver(a0, u0, p0, dt, cmk, N, L, T, rho));
            shared_ptr<TubeFlowSolidSolver> solid(new TubeFlowSolidSolver(a0, cmk, p0, rho, L, N));

            // Convergence measures
            std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
            convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >(new std::list<std::shared_ptr<ConvergenceMeasure> > );

            // Add convergence measure
            convergenceMeasures->push_back(std::shared_ptr<ConvergenceMeasure> (new RelativeConvergenceMeasure(0, false, tol)));

            // Create the fsi solver
            fsi = new FsiSolver(fluid, solid, convergenceMeasures, parallel, extrapolationOrder);
        }

        virtual void TearDown() {
            delete fsi;
        }

        FsiSolver * fsi;
};

INSTANTIATE_TEST_CASE_P(testParameters, FsiSolverTest, ::testing::Bool());

TEST_P(FsiSolverTest, object)
{
    ASSERT_TRUE(true);
}

TEST_P(FsiSolverTest, initTimeStep)
{
    ASSERT_FALSE(fsi->init);
    ASSERT_FALSE(fsi->fluid->init);
    ASSERT_FALSE(fsi->solid->init);

    fsi->initTimeStep();

    ASSERT_TRUE(fsi->init);
    ASSERT_TRUE(fsi->fluid->init);
    ASSERT_TRUE(fsi->solid->init);
}

TEST_P(FsiSolverTest, evaluate)
{
    bool parallel = GetParam();

    fsi::vector input(5), output(5), R(5);

    if (not parallel)
        input = fsi->solid->data.col(0);

    if (parallel) {
        input.resize(10);
        output.resize(10);
        R.resize(10);
        input.head(5) = fsi->solid->data.col(0);
        input.tail(5) = fsi->fluid->data.col(0);
    }

    output.setZero();
    R.setZero();

    fsi->initTimeStep();

    fsi->evaluate(input, output, R);
}
