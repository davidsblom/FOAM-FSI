
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <unsupported/Eigen/NumericalDiff>

#include "AitkenPostProcessing.H"
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

class AitkenPostProcessingTest : public::testing::Test
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
        int N = 10;
        double L = 1;
        double T = 10;
        double dx = L / N;
        double rho = 1.225;
        double E = 490;
        double h = 1.0e-3;
        double cmk = std::sqrt( E * h / (2 * rho * r0) );
        double c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
        double kappa = c0 / u0;
        double tau = u0 * dt / L;

        // Computational settings
        double tol = 1.0e-7;
        int maxIter = 200;
        double initialRelaxation = 1.0e-3;
        int reuseInformationStartingFromTimeIndex = 0;
        int minIter = 1;
        int extrapolation = 2;

        // Parametrized settings
        bool parallel = false;
        int nbReuse = 0;
        int maxUsedIterations = N;

        if ( parallel )
            maxUsedIterations *= 2;

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        shared_ptr<TubeFlowFluidSolver> fluid( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
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

        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, minIter ) ) );
        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 1, tol ) ) );

        shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

        shared_ptr<PostProcessing> postProcessing;
        postProcessing = shared_ptr<PostProcessing>( new AitkenPostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

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

TEST_F( AitkenPostProcessingTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( AitkenPostProcessingTest, initTimeStep )
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

TEST_F( AitkenPostProcessingTest, iqn_evaluate_residual )
{
    solver->initTimeStep();

    int N = 10;

    fsi::vector input( N ), output( N ), R( N );
    input.head( 10 ) = solver->fsi->solid->data.col( 0 );
    output.setZero();
    R.setZero();

    solver->postProcessing->fsi->evaluate( input, output, R );
}

TEST_F( AitkenPostProcessingTest, timestep )
{
    solver->solveTimeStep();

    ASSERT_TRUE( solver->fsi->allConverged );
    ASSERT_TRUE( solver->fsi->fluid->isRunning() );
}

TEST_F( AitkenPostProcessingTest, monolithic )
{
    for ( int i = 0; i < 100; i++ )
    {
        solver->solveTimeStep();
        monolithicSolver->solveTimeStep();

        double tol = 1.0e-5;

        if ( i < 99 )
        {
            ASSERT_TRUE( solver->fsi->fluid->isRunning() );
            ASSERT_TRUE( solver->fsi->solid->isRunning() );
        }
        else
            ASSERT_FALSE( solver->fsi->fluid->isRunning() );

        ASSERT_TRUE( solver->fsi->allConverged );
        ASSERT_NEAR( solver->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
        ASSERT_NEAR( solver->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
        ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
        ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );
    }
}
