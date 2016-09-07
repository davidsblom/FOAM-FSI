
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AndersonPostProcessing.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
#include "MultiLevelSpaceMappingSolver.H"
#include "OutputSpaceMapping.H"
#include "RelativeConvergenceMeasure.H"
#include "SpaceMappingSolver.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class MultiLevelOutputSpaceMappingSolverParametrizedTest : public TestWithParam< std::tr1::tuple<int, bool> >
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
        int N = 10;
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
        int maxIter = 200;
        scalar initialRelaxation = 1.0e-3;
        bool scaling = false;
        scalar singularityLimit = 1.0e-13;
        int reuseInformationStartingFromTimeIndex = 0;
        int order = 1;
        scalar beta = 1;
        bool updateJacobian = false;
        bool parallel = false;
        int extrapolation = 2;
        int minIter = 2;
        int nbReuse = 2;

        // Parametrized settings
        int couplingGridSize = std::tr1::get<0>( GetParam() );
        bool synchronization = std::tr1::get<1>( GetParam() );

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        // Initialize variables
        shared_ptr<TubeFlowFluidSolver> fluid;
        shared_ptr<TubeFlowSolidSolver> solid;
        shared_ptr<TubeFlowFluidSolver> fineModelFluid;
        shared_ptr<TubeFlowSolidSolver> fineModelSolid;
        shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        shared_ptr<AndersonPostProcessing> postProcessing;
        shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > solvers;
        shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

        solvers = shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > ( new std::deque<shared_ptr<SpaceMappingSolver> >() );
        models = shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > ( new std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );

        fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
        fineModelSolid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

        // Coarse model

        int maxUsedIterations = N;

        if ( parallel )
            maxUsedIterations *= 2;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
        shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
        shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
        shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 1.0e-3 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 1.0e-3 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        fluid.reset();
        solid.reset();
        multiLevelFluidSolver.reset();
        multiLevelSolidSolver.reset();
        convergenceMeasures.reset();
        multiLevelFsiSolver.reset();
        postProcessing.reset();

        // Fine model

        maxUsedIterations = 2 * N;

        if ( parallel )
            maxUsedIterations *= 2;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, 2 * N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, 2 * N ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 1 ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 1 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 1.0e-2 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 1.0e-2 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        shared_ptr<OutputSpaceMapping> outputSpaceMapping( new OutputSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

        // Create manifold mapping solver
        shared_ptr<SpaceMappingSolver > outputSpaceMappingSolver( new SpaceMappingSolver( fineModel, coarseModel, outputSpaceMapping ) );

        solvers->push_back( outputSpaceMappingSolver );
        models->push_back( coarseModel );
        models->push_back( fineModel );

        shared_ptr<SpaceMappingSolver> coarseModelSolver = outputSpaceMappingSolver;

        fluid.reset();
        solid.reset();
        multiLevelFluidSolver.reset();
        multiLevelSolidSolver.reset();
        convergenceMeasures.reset();
        multiLevelFsiSolver.reset();
        postProcessing.reset();

        // Fine model 2

        maxUsedIterations = couplingGridSize;

        if ( parallel )
            maxUsedIterations *= 2;

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 2 ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 2 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, 1 ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        outputSpaceMapping = shared_ptr<OutputSpaceMapping> ( new OutputSpaceMapping( fineModel, coarseModelSolver, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

        // Create manifold mapping solver
        outputSpaceMappingSolver = shared_ptr<SpaceMappingSolver> ( new SpaceMappingSolver( fineModel, coarseModelSolver, outputSpaceMapping ) );

        solvers->push_back( outputSpaceMappingSolver );
        models->push_back( fineModel );

        // Create multi level manifold mapping solver

        solver = new MultiLevelSpaceMappingSolver( solvers, models, synchronization );

        // Monolithic solver
        monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho );
    }

    virtual void TearDown()
    {
        delete solver;
        delete monolithicSolver;
    }

    MonolithicFsiSolver * monolithicSolver;
    MultiLevelSpaceMappingSolver * solver;
};

INSTANTIATE_TEST_CASE_P( testParameters, MultiLevelOutputSpaceMappingSolverParametrizedTest, ::testing::Combine( Values( 20, 40 ), Bool() ) );

TEST_P( MultiLevelOutputSpaceMappingSolverParametrizedTest, object )
{
    ASSERT_TRUE( true );
}


TEST_P( MultiLevelOutputSpaceMappingSolverParametrizedTest, run )
{
    solver->solveTimeStep();

    ASSERT_TRUE( solver->models->at( solver->models->size() - 1 )->fsi->allConverged );
    ASSERT_TRUE( solver->models->at( solver->models->size() - 1 )->fsi->fluid->isRunning() );
    ASSERT_EQ( solver->solvers->at( 0 )->spaceMapping->timeIndex, 1 );
    ASSERT_EQ( solver->solvers->at( 1 )->spaceMapping->timeIndex, 1 );
}


TEST_P( MultiLevelOutputSpaceMappingSolverParametrizedTest, monolithic )
{
    bool synchronization = std::tr1::get<1>( GetParam() );

    for ( int i = 0; i < 5; i++ )
    {
        solver->solveTimeStep();
        monolithicSolver->solveTimeStep();

        scalar tol = 1.0e-5;

        if ( i < 99 )
            ASSERT_TRUE( solver->models->at( solver->models->size() - 1 )->fsi->fluid->isRunning() );
        else
            ASSERT_FALSE( solver->models->at( solver->models->size() - 1 )->fsi->fluid->isRunning() );

        ASSERT_TRUE( solver->models->at( 2 )->fsi->allConverged );
        ASSERT_NEAR( solver->models->at( solver->models->size() - 1 )->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
        ASSERT_NEAR( solver->models->at( solver->models->size() - 1 )->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
        ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
        ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );

        // Verify that the coarse models are synchronized with the fine model
        if ( synchronization )
        {
            for ( auto && model : *(solver->models) )
                ASSERT_NEAR( model->fsi->x.norm(), solver->solvers->at( solver->solvers->size() - 1 )->fsi->x.norm(), 1.0e-13 );
        }
    }
}
