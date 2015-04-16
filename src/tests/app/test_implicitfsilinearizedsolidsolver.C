
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AndersonPostProcessing.H"
#include "ConvergenceMeasure.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowLinearizedFluidSolver.H"
#include "TubeFlowLinearizedSolidSolver.H"
#include "gtest/gtest.h"
#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class ImplicitFsiLinearizedSolidSolverParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, bool> >
{
protected:

    virtual void SetUp()
    {
        // Physical settings
        double r0 = 0.2;
        double a0 = M_PI * r0 * r0;
        double u0 = 0.1;
        double p0 = 0;
        double dt = 0.1;
        int N = 5;
        double L = 1;
        double T = 10;
        double dx = L / N;
        double rho = 1.225;
        double E = 490;
        double G = E;
        double nu = 0.5;
        double h = 1.0e-3;
        double cmk = std::sqrt( E * h / (2 * rho * r0) );
        double c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
        double kappa = c0 / u0;
        double tau = u0 * dt / L;

        // Computational settings
        double tol = 1.0e-7;
        int maxIter = 50;
        double initialRelaxation = 1.0e-3;
        double singularityLimit = 1.0e-11;
        int reuseInformationStartingFromTimeIndex = 0;
        bool scaling = false;
        double beta = 1;
        bool updateJacobian = false;
        int extrapolation = 2;

        // Parametrized settings
        bool parallel = std::tr1::get<0>( GetParam() );
        int nbReuse = std::tr1::get<1>( GetParam() );
        bool convergenceMeasureTraction = std::tr1::get<2>( GetParam() );

        int maxUsedIterations = N;

        if ( parallel )
            maxUsedIterations *= 2;

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        shared_ptr<TubeFlowLinearizedFluidSolver> fluid( new TubeFlowLinearizedFluidSolver( N, p0, r0, u0, rho, E, h, T, dt, L ) );
        shared_ptr<TubeFlowLinearizedSolidSolver> solid( new TubeFlowLinearizedSolidSolver( N, nu, rho, h, L, dt, G, E, r0 ) );

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

        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, tol ) ) );

        if ( parallel || convergenceMeasureTraction )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 1, tol ) ) );

        shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );
        shared_ptr<AndersonPostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );
        solver = new ImplicitMultiLevelFsiSolver( fsi, postProcessing );
        monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, N, L, T, rho );
    }

    virtual void TearDown()
    {
        delete solver;
        delete monolithicSolver;
    }

    ImplicitMultiLevelFsiSolver * solver;
    MonolithicFsiSolver * monolithicSolver;
};

INSTANTIATE_TEST_CASE_P( testParameters, ImplicitFsiLinearizedSolidSolverParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Bool() ) );

TEST_P( ImplicitFsiLinearizedSolidSolverParametrizedTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( ImplicitFsiLinearizedSolidSolverParametrizedTest, initTimeStep )
{
    ASSERT_FALSE( solver->init );
    ASSERT_FALSE( solver->fsi->init );
    ASSERT_FALSE( solver->fsi->fluid->init );
    ASSERT_FALSE( solver->fsi->solid->init );
    ASSERT_FALSE( solver->postProcessing->fsi->fluid->init );
    ASSERT_FALSE( solver->postProcessing->fsi->init );
    ASSERT_FALSE( solver->postProcessing->fsi->fluid->init );

    solver->initTimeStep();

    ASSERT_TRUE( solver->init );
    ASSERT_TRUE( solver->fsi->init );
    ASSERT_TRUE( solver->fsi->fluid->init );
    ASSERT_TRUE( solver->fsi->solid->init );
    ASSERT_TRUE( solver->postProcessing->fsi->fluid->init );
    ASSERT_TRUE( solver->postProcessing->fsi->fluid->init );
}

TEST_P( ImplicitFsiLinearizedSolidSolverParametrizedTest, residual )
{
    bool parallel = std::tr1::get<0>( GetParam() );
    int N = 5;

    if ( parallel )
        N = 10;

    fsi::vector input( N ), output( N ), R( N );
    input.head( 5 ) = solver->fsi->solid->data.col( 0 );
    output.setZero();
    R.setZero();

    solver->initTimeStep();

    solver->fsi->evaluate( input, output, R );
}

TEST_P( ImplicitFsiLinearizedSolidSolverParametrizedTest, iqn_evaluate_residual )
{
    solver->initTimeStep();

    bool parallel = std::tr1::get<0>( GetParam() );
    int N = 5;

    if ( parallel )
        N = 10;

    fsi::vector input( N ), output( N ), R( N );
    input.head( 5 ) = solver->fsi->solid->data.col( 0 );
    output.setZero();
    R.setZero();

    solver->postProcessing->fsi->evaluate( input, output, R );
}

TEST_P( ImplicitFsiLinearizedSolidSolverParametrizedTest, timestep )
{
    solver->solveTimeStep();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_TRUE( solver->fsi->fluid->isRunning() );
}

TEST_P( ImplicitFsiLinearizedSolidSolverParametrizedTest, run )
{
    solver->run();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_FALSE( solver->fsi->fluid->isRunning() );
}
