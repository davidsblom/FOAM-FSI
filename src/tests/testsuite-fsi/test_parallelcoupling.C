
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ConvergenceMeasure.H"
#include "FsiSolver.H"
#include "AndersonPostProcessing.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;
using ::testing::TestWithParam;
using ::testing::Bool;
using ::testing::Values;
using ::testing::Combine;

class parallelCouplingParametrizedTest : public TestWithParam<int>
{
protected:

    virtual void SetUp()
    {
        // Physical settings
        scalar r0 = 0.2;
        scalar a0 = M_PI * r0 * r0;
        scalar u0 = 0.1;
        scalar p0 = 0;
        scalar dt = 0.1;
        int N = GetParam();
        scalar L = 1;
        scalar T = 10;
        scalar dx = L / N;
        scalar rho = 1.225;
        scalar E = 490;
        scalar h = 1.0e-3;
        scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
        scalar c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
        scalar kappa = c0 / u0;
        scalar tau = u0 * dt / L;

        // Computational settings
        scalar tol = 1.0e-5;
        int maxIter = 50;
        scalar initialRelaxation = 1.0e-3;
        scalar singularityLimit = 1.0e-11;
        int reuseInformationStartingFromTimeIndex = 10;
        int maxUsedIterations = 2 * N;
        bool scaling = true;
        scalar beta = 1;
        bool updateJacobian = false;

        // Parametrized settings
        bool parallel = true;
        int nbReuse = 0;
        int extrapolation = 0;
        int minIter = 1;

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        shared_ptr<TubeFlowSolidSolver> solid( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        shared_ptr<RBFFunctionInterface> rbfFunction;
        shared_ptr<RBFInterpolation> rbfInterpolator;
        shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
        shared_ptr<RBFCoarsening> rbfInterpToMesh;

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        // Convergence measures
        std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 1, false, tol ) ) );

        shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );
        shared_ptr<AndersonPostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );
        solver = new ImplicitMultiLevelFsiSolver( fsi, postProcessing );
    }

    virtual void TearDown()
    {
        delete solver;
    }

    ImplicitMultiLevelFsiSolver * solver;
    std::shared_ptr<TubeFlowFluidSolver> fluid;
};

INSTANTIATE_TEST_CASE_P( parallelCoupling, parallelCouplingParametrizedTest, Values( 5, 100 ) );

TEST_P( parallelCouplingParametrizedTest, timeStep )
{
    int N = GetParam();

    solver->solveTimeStep();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_TRUE( solver->fsi->fluid->isRunning() );

    if ( N == 5 )
    {
        ASSERT_LE( solver->fsi->nbIter, 12 );
        ASSERT_LE( fluid->nbRes, 60 );
        ASSERT_LE( fluid->nbJac, 36 );
    }

    if ( N == 100 )
    {
        ASSERT_LE( solver->fsi->nbIter, 20 );
        ASSERT_LE( fluid->nbRes, 100 );
        ASSERT_LE( fluid->nbJac, 60 );
    }
}

TEST_P( parallelCouplingParametrizedTest, timeSteps )
{
    int N = GetParam();

    solver->solveTimeStep();

    if ( N == 100 )
    {
        ASSERT_LE( fluid->nbRes, 100 );
        ASSERT_LE( fluid->nbJac, 60 );
    }

    solver->solveTimeStep();

    if ( N == 100 )
    {
        ASSERT_LE( fluid->nbRes, 215 );
        ASSERT_LE( fluid->nbJac, 129 );
    }

    solver->solveTimeStep();

    if ( N == 100 )
    {
        ASSERT_LE( fluid->nbRes, 335 );
        ASSERT_LE( fluid->nbJac, 201 );
    }

    solver->solveTimeStep();

    if ( N == 100 )
    {
        ASSERT_LE( fluid->nbRes, 474 );
        ASSERT_LE( fluid->nbJac, 285 );
    }

    solver->solveTimeStep();

    if ( N == 100 )
    {
        ASSERT_LE( fluid->nbRes, 593 );
        ASSERT_LE( fluid->nbJac, 357 );
    }

    solver->solveTimeStep();
}

TEST_P( parallelCouplingParametrizedTest, run )
{
    solver->run();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_FALSE( solver->fsi->fluid->isRunning() );
}
