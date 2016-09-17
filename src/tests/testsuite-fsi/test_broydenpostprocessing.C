
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <unsupported/Eigen/NumericalDiff>

#include "BroydenPostProcessing.H"
#include "ConvergenceMeasure.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using ::testing::TestWithParam;
using ::testing::Bool;
using ::testing::Values;
using ::testing::Combine;

class BroydenPostProcessingParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int> >{
    protected:
        virtual void SetUp() {
            // Physical settings
            scalar r0 = 0.2;
            scalar a0 = M_PI * r0 * r0;
            scalar u0 = 0.1;
            scalar p0 = 0;
            scalar dt = 0.1;
            int N = 10;
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
            int maxIter = 200;
            scalar initialRelaxation = 1.0e-3;
            int reuseInformationStartingFromTimeIndex = 0;
            int minIter = 1;
            int extrapolation = 2;
            scalar singularityLimit = 1.0e-11;

            // Parametrized settings
            bool parallel = std::tr1::get<0>(GetParam());
            int nbReuse = std::tr1::get<1>(GetParam());
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

            shared_ptr<PostProcessing> postProcessing;
            postProcessing = shared_ptr<PostProcessing>(new BroydenPostProcessing(fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex));

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

INSTANTIATE_TEST_CASE_P(testParameters, BroydenPostProcessingParametrizedTest, ::testing::Combine(Bool(), Values(0, 1, 4, 8)));

TEST_P(BroydenPostProcessingParametrizedTest, object)
{
    ASSERT_TRUE(true);
}

TEST_P(BroydenPostProcessingParametrizedTest, initTimeStep)
{
    ASSERT_FALSE(solver->init);
    ASSERT_FALSE(solver->fsi->init);
    ASSERT_FALSE(solver->fsi->fluid->init);
    ASSERT_FALSE(solver->fsi->solid->init);
    ASSERT_FALSE(solver->postProcessing->fsi->fluid->init);
    ASSERT_FALSE(solver->postProcessing->fsi->init);
    ASSERT_FALSE(solver->postProcessing->fsi->fluid->init);

    solver->initTimeStep();

    ASSERT_TRUE(solver->init);
    ASSERT_TRUE(solver->fsi->init);
    ASSERT_TRUE(solver->fsi->fluid->init);
    ASSERT_TRUE(solver->fsi->solid->init);
    ASSERT_TRUE(solver->postProcessing->fsi->fluid->init);
    ASSERT_TRUE(solver->postProcessing->fsi->fluid->init);
}

TEST_P(BroydenPostProcessingParametrizedTest, iqn_evaluate_residual)
{
    solver->initTimeStep();

    bool parallel = std::tr1::get<0>(GetParam());
    int N = 10;

    if (parallel)
        N = 20;

    fsi::vector input(N), output(N), R(N);
    input.head(10) = solver->fsi->solid->data.col(0);
    output.setZero();
    R.setZero();

    solver->postProcessing->fsi->evaluate(input, output, R);
}

TEST_P(BroydenPostProcessingParametrizedTest, timestep)
{
    solver->solveTimeStep();

    ASSERT_TRUE(solver->fsi->allConverged);
    ASSERT_TRUE(solver->fsi->fluid->isRunning());
}

TEST_P(BroydenPostProcessingParametrizedTest, monolithic)
{
    for (int i = 0; i < 100; i++) {
        solver->solveTimeStep();
        monolithicSolver->solveTimeStep();

        scalar tol = 1.0e-5;

        if (i < 99) {
            ASSERT_TRUE(solver->fsi->fluid->isRunning());
            ASSERT_TRUE(solver->fsi->solid->isRunning());
        }else
            ASSERT_FALSE(solver->fsi->fluid->isRunning());

        ASSERT_TRUE(solver->fsi->allConverged);
        ASSERT_NEAR(solver->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol);
        ASSERT_NEAR(solver->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol);
        ASSERT_TRUE(monolithicSolver->an.norm() > 0);
        ASSERT_TRUE(monolithicSolver->pn.norm() > 0);
    }
}

TEST_P(BroydenPostProcessingParametrizedTest, reuse)
{
    solver->run();

    bool parallel = std::tr1::get<0>(GetParam());
    int nbReuse = std::tr1::get<1>(GetParam());

    if (!parallel && nbReuse == 0)
        ASSERT_LE(solver->fsi->nbIter, 1072);

    if (!parallel && nbReuse == 1)
        ASSERT_LE(solver->fsi->nbIter, 703);

    if (!parallel && nbReuse == 4)
        ASSERT_LE(solver->fsi->nbIter, 582);

    if (!parallel && nbReuse == 8)
        ASSERT_LE(solver->fsi->nbIter, 530);

    if (parallel && nbReuse == 0)
        ASSERT_LE(solver->fsi->nbIter, 2069);

    if (parallel && nbReuse == 1)
        ASSERT_LE(solver->fsi->nbIter, 1480);

    if (parallel && nbReuse == 4)
        ASSERT_LE(solver->fsi->nbIter, 1124);

    if (parallel && nbReuse == 8)
        ASSERT_LE(solver->fsi->nbIter, 993);
}

TEST_P(BroydenPostProcessingParametrizedTest, numberOfColumnsVIQN)
{
    int nbReuse = std::tr1::get<1>(GetParam());
    int nbIterFirstTimeStep = 0;

    for (int i = 0; i < 3; i++) {
        solver->initTimeStep();
        solver->solve();

        // Determine the number of columns of the V and W matrices
        int nbResiduals = solver->postProcessing->residuals.size();

        // Include information from previous optimization solves
        for (auto sols : solver->postProcessing->solsList)
            nbResiduals += sols.size();

        // Include information from previous stages
        for (auto solsList : solver->postProcessing->solsStageList)
            for (auto sols : solsList)
                nbResiduals += sols.size();

        // Include information from previous time steps
        for (auto solsStageList : solver->postProcessing->solsTimeList)
            for (auto solsList : solsStageList)
                for (auto sols : solsList)
                    nbResiduals += sols.size();

        if (i == 0)
            nbIterFirstTimeStep = solver->fsi->iter;

        if (nbReuse == 0)
            ASSERT_EQ(nbResiduals, solver->fsi->iter);

        if (nbReuse == 1 && i < 2)
            ASSERT_EQ(nbResiduals, solver->fsi->nbIter);

        if (nbReuse == 1 && i == 2)
            ASSERT_EQ(nbResiduals, solver->fsi->nbIter - nbIterFirstTimeStep);

        if (nbReuse > 1)
            ASSERT_EQ(nbResiduals, solver->fsi->nbIter);

        solver->finalizeTimeStep();
    }
}
