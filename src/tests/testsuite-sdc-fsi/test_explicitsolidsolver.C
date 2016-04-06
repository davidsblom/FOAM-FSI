
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AndersonPostProcessing.H"
#include "ConvergenceMeasure.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MultiLevelFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowLinearizedFluidSolver.H"
#include "TubeFlowExplicitLinearSolidSolver.H"
#include "gtest/gtest.h"
#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;
using ::testing::TestWithParam;
using ::testing::Bool;
using ::testing::Values;
using ::testing::Combine;

class ExplicitLinearSolidSolverTest : public TestWithParam< std::tr1::tuple<bool, int, bool> >
{
    protected:
        virtual void SetUp()
        {
            // Physical settings
            scalar r0 = 0.2;
            scalar u0 = 0.1;
            scalar p0 = 0;
            scalar dt = 0.01;
            int N = 5;
            scalar L = 1;
            scalar T = 10;
            scalar dx = L / N;
            scalar rho = 1.225;
            scalar E = 490;
            scalar G = E;
            scalar nu = 0.5;
            scalar h = 1.0e-3;
            scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
            scalar c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
            scalar kappa = c0 / u0;

            // Computational settings
            scalar tol = 1.0e-7;
            int maxIter = 50;
            scalar initialRelaxation = 1.0e-3;
            scalar singularityLimit = 1.0e-11;
            int reuseInformationStartingFromTimeIndex = 0;
            bool scaling = false;
            scalar beta = 1;
            bool updateJacobian = false;
            int extrapolation = 2;

            // Parametrized settings
            bool parallel = std::tr1::get<0>( GetParam() );
            int nbReuse = std::tr1::get<1>( GetParam() );
            bool convergenceMeasureTraction = std::tr1::get<2>( GetParam() );

            int maxUsedIterations = N;

            if ( parallel )
                maxUsedIterations *= 2;

            ASSERT_NEAR( kappa, 10, 1.0e-13 );
            ASSERT_TRUE( dx > 0 );

            shared_ptr<TubeFlowLinearizedFluidSolver> fluid( new TubeFlowLinearizedFluidSolver( N, p0, r0, u0, rho, E, h, T, dt, L ) );
            shared_ptr<TubeFlowExplicitLinearSolidSolver> solid( new TubeFlowExplicitLinearSolidSolver( N, nu, rho, h, L, dt, G, E, r0, p0 ) );

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

            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, false, tol ) ) );

            if ( parallel || convergenceMeasureTraction )
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
};

INSTANTIATE_TEST_CASE_P( testParameters, ExplicitLinearSolidSolverTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Bool() ) );

TEST_P( ExplicitLinearSolidSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( ExplicitLinearSolidSolverTest, initTimeStep )
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

TEST_P( ExplicitLinearSolidSolverTest, residual )
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

TEST_P( ExplicitLinearSolidSolverTest, iqn_evaluate_residual )
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

TEST_P( ExplicitLinearSolidSolverTest, timestep )
{
    solver->solveTimeStep();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_TRUE( solver->fsi->fluid->isRunning() );
}

TEST_P( ExplicitLinearSolidSolverTest, run )
{
    solver->run();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_FALSE( solver->fsi->fluid->isRunning() );
}
