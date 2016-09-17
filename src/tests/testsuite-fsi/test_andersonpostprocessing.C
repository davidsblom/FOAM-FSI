
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ConvergenceMeasure.H"
#include "FsiSolver.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"
#include "AndersonPostProcessing.H"
#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;
using ::testing::TestWithParam;
using ::testing::Bool;
using ::testing::Values;
using ::testing::Combine;

class AndersonPostProcessingParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, bool> >{
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
            scalar dx = L / N;
            scalar rho = 1.225;
            scalar E = 490;
            scalar h = 1.0e-3;
            scalar cmk = std::sqrt(E * h / (2 * rho * r0));
            scalar c0 = std::sqrt(cmk * cmk - p0 / (2 * rho));
            scalar kappa = c0 / u0;
            scalar tau = u0 * dt / L;

            // Computational settings
            scalar tol = 1.0e-7;
            int maxIter = 50;
            scalar initialRelaxation = 1.0e-3;
            scalar singularityLimit = 1.0e-11;
            int reuseInformationStartingFromTimeIndex = 0;
            bool scaling = false;
            scalar beta = 1;

            // Parametrized settings
            bool parallel = std::tr1::get<0>(GetParam());
            int nbReuse = std::tr1::get<1>(GetParam());
            int extrapolation = std::tr1::get<2>(GetParam());
            int minIter = std::tr1::get<3>(GetParam());
            bool updateJacobian = std::tr1::get<4>(GetParam());

            int maxUsedIterations = N;

            if (parallel)
                maxUsedIterations *= 2;

            ASSERT_NEAR(tau, 0.01, 1.0e-13);
            ASSERT_NEAR(kappa, 10, 1.0e-13);
            ASSERT_TRUE(dx > 0);

            shared_ptr<TubeFlowFluidSolver> fluid(new TubeFlowFluidSolver(a0, u0, p0, dt, cmk, N, L, T, rho));
            shared_ptr<TubeFlowSolidSolver> solid(new TubeFlowSolidSolver(a0, cmk, p0, rho, L, N));

            shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
            shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
            shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
            shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

            rbfFunction = shared_ptr<rbf::RBFFunctionInterface>(new rbf::TPSFunction());
            rbfInterpolator = shared_ptr<rbf::RBFInterpolation>(new rbf::RBFInterpolation(rbfFunction));
            rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> (new rbf::RBFCoarsening(rbfInterpolator));

            rbfFunction = shared_ptr<rbf::RBFFunctionInterface>(new rbf::TPSFunction());
            rbfInterpolator = shared_ptr<rbf::RBFInterpolation>(new rbf::RBFInterpolation(rbfFunction));
            rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> (new rbf::RBFCoarsening(rbfInterpolator));

            shared_ptr<MultiLevelSolver> fluidSolver(new MultiLevelSolver(fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0));

            rbfFunction = shared_ptr<rbf::RBFFunctionInterface>(new rbf::TPSFunction());
            rbfInterpolator = shared_ptr<rbf::RBFInterpolation>(new rbf::RBFInterpolation(rbfFunction));
            rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> (new rbf::RBFCoarsening(rbfInterpolator));

            rbfFunction = shared_ptr<rbf::RBFFunctionInterface>(new rbf::TPSFunction());
            rbfInterpolator = shared_ptr<rbf::RBFInterpolation>(new rbf::RBFInterpolation(rbfFunction));
            rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> (new rbf::RBFCoarsening(rbfInterpolator));

            shared_ptr<MultiLevelSolver> solidSolver(new MultiLevelSolver(solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0));

            // Convergence measures
            std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
            convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >(new std::list<std::shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back(std::shared_ptr<ConvergenceMeasure>(new MinIterationConvergenceMeasure(0, false, minIter)));
            convergenceMeasures->push_back(std::shared_ptr<ConvergenceMeasure>(new RelativeConvergenceMeasure(0, false, tol)));

            if (parallel)
                convergenceMeasures->push_back(std::shared_ptr<ConvergenceMeasure>(new RelativeConvergenceMeasure(1, false, tol)));

            shared_ptr<MultiLevelFsiSolver> fsi(new MultiLevelFsiSolver(fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation));
            shared_ptr<AndersonPostProcessing> postProcessing(new AndersonPostProcessing(fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian));
            solver = new ImplicitMultiLevelFsiSolver(fsi, postProcessing);
            monolithicSolver = new MonolithicFsiSolver(a0, u0, p0, dt, cmk, N, L, T, rho);
        }

        virtual void TearDown() {
            delete solver;
            delete monolithicSolver;
        }

        ImplicitMultiLevelFsiSolver * solver;
        MonolithicFsiSolver * monolithicSolver;
};

INSTANTIATE_TEST_CASE_P(testParameters, AndersonPostProcessingParametrizedTest, ::testing::Combine(Bool(), Values(0, 1, 4), Values(0, 1, 2), Values(1), Bool()));

TEST_P(AndersonPostProcessingParametrizedTest, object)
{
    ASSERT_TRUE(true);
}

TEST_P(AndersonPostProcessingParametrizedTest, run)
{
    solver->run();
    monolithicSolver->run();

    scalar tol = 1.0e-5;
    ASSERT_TRUE(solver->fsi->allConverged);
    ASSERT_NEAR(solver->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol);
    ASSERT_NEAR(solver->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol);
    ASSERT_TRUE(monolithicSolver->an.norm() > 0);
    ASSERT_TRUE(monolithicSolver->pn.norm() > 0);
}
