
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
#include "OutputSpaceMapping.H"
#include "RelativeConvergenceMeasure.H"
#include "SpaceMappingSolver.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "MultiLevelFsiSolver.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class OutputSpaceMappingSolverParametrizedTest : public TestWithParam< std::tr1::tuple<int, int, int, int, int> >
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
        int maxIter = 50;
        scalar initialRelaxation = 1.0e-3;
        int reuseInformationStartingFromTimeIndex = 0;
        scalar singularityLimit = 1.0e-13;
        bool scaling = false;
        scalar beta = 1;
        bool updateJacobian = false;

        // Parametrized settings
        bool parallel = false;
        int nbReuse = std::tr1::get<0>( GetParam() );
        int extrapolation = std::tr1::get<1>( GetParam() );
        int minIter = std::tr1::get<2>( GetParam() );
        int couplingGridSize = std::tr1::get<3>( GetParam() );
        int order = std::tr1::get<4>( GetParam() );

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        // Initialize variables
        shared_ptr<TubeFlowFluidSolver> fluid;
        shared_ptr<TubeFlowSolidSolver> solid;
        shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        shared_ptr<AndersonPostProcessing> postProcessing;
        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

        int maxUsedIterations = couplingGridSize;

        if ( parallel )
            maxUsedIterations *= 2;

        // Fine model

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

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

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 1 ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 1 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, 1 ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Coarse model

        maxUsedIterations = N;

        if ( parallel )
            maxUsedIterations *= 2;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->solid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 0.1 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 0.1 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        shared_ptr<OutputSpaceMapping> outputSpaceMapping( new OutputSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

        // Create manifold mapping solver
        solver = new SpaceMappingSolver( fineModel, coarseModel, outputSpaceMapping );

        // Monolithic solver
        monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho );

        assert( fluid );
        assert( solid );
        assert( multiLevelFluidSolver );
        assert( multiLevelSolidSolver );
        assert( convergenceMeasures );
        assert( multiLevelFsiSolver );
        assert( postProcessing );
    }

    virtual void TearDown()
    {
        delete solver;
        delete monolithicSolver;
    }

    MonolithicFsiSolver * monolithicSolver;
    SpaceMappingSolver * solver;
};

INSTANTIATE_TEST_CASE_P( testParameters, OutputSpaceMappingSolverParametrizedTest, ::testing::Combine( Values( 0, 1, 4 ), Values( 2 ), Values( 3 ), Values( 10, 20 ), Values( 0, 1, 2 ) ) );

TEST_P( OutputSpaceMappingSolverParametrizedTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( OutputSpaceMappingSolverParametrizedTest, run )
{
    for ( int i = 0; i < 10; i++ )
    {
        solver->solveTimeStep();
        monolithicSolver->solveTimeStep();
    }

    scalar tol = 1.0e-5;
    ASSERT_TRUE( solver->fineModel->fsi->allConverged );
    ASSERT_TRUE( solver->coarseModel->fsi->allConverged );
    ASSERT_NEAR( solver->fineModel->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
    ASSERT_NEAR( solver->fineModel->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
    ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
    ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );
}

TEST_P( OutputSpaceMappingSolverParametrizedTest, iterations )
{
    solver->run();

    int nbReuse = std::tr1::get<0>( GetParam() );
    int extrapolation = std::tr1::get<1>( GetParam() );
    int minIter = std::tr1::get<2>( GetParam() );
    int couplingGridSize = std::tr1::get<3>( GetParam() );
    int order = std::tr1::get<4>( GetParam() );

    if ( couplingGridSize == 20 && nbReuse == 0 && extrapolation == 2 && minIter == 3 && order == 0 )
    {
        ASSERT_LE( solver->fineModel->fsi->nbIter, 874 );
        ASSERT_LE( solver->coarseModel->fsi->nbIter, 4049 );
    }

    if ( couplingGridSize == 20 && nbReuse == 1 && extrapolation == 2 && minIter == 3 && order == 0 )
    {
        ASSERT_LE( solver->fineModel->fsi->nbIter, 874 );
        ASSERT_LE( solver->coarseModel->fsi->nbIter, 3713 );
    }

    if ( couplingGridSize == 20 && nbReuse == 4 && extrapolation == 2 && minIter == 3 && order == 0 )
    {
        ASSERT_LE( solver->fineModel->fsi->nbIter, 874 );
        ASSERT_LE( solver->coarseModel->fsi->nbIter, 3713 );
    }

    if ( couplingGridSize == 20 && nbReuse == 0 && extrapolation == 2 && minIter == 3 && order == 1 )
    {
        ASSERT_LE( solver->fineModel->fsi->nbIter, 710 );
        ASSERT_LE( solver->coarseModel->fsi->nbIter, 5198 );
    }

    if ( couplingGridSize == 20 && nbReuse == 1 && extrapolation == 2 && minIter == 3 && order == 1 )
    {
        ASSERT_LE( solver->fineModel->fsi->nbIter, 668 );
        ASSERT_LE( solver->coarseModel->fsi->nbIter, 5252 );
    }

    if ( couplingGridSize == 20 && nbReuse == 4 && extrapolation == 2 && minIter == 3 && order == 1 )
    {
        ASSERT_LE( solver->fineModel->fsi->nbIter, 574 );
        ASSERT_LE( solver->coarseModel->fsi->nbIter, 4132 );
    }
}

class OutputSpaceMappingSolverTest : public::testing::Test
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
        int maxIter = 50;
        scalar initialRelaxation = 1.0e-3;
        int reuseInformationStartingFromTimeIndex = 0;
        scalar singularityLimit = 1.0e-11;
        bool scaling = false;
        scalar beta = 1;
        bool updateJacobian = false;

        // Parametrized settings
        bool parallel = false;
        int nbReuse = 0;
        int extrapolation = 0;
        int minIter = 3;
        int couplingGridSize = 100;
        int order = 1;

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        // Initialize variables
        shared_ptr<TubeFlowSolidSolver> solid;
        shared_ptr<TubeFlowFluidSolver> fluid;
        shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        shared_ptr<AndersonPostProcessing> postProcessing;
        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

        // Fine model

        int maxUsedIterations = couplingGridSize;

        fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

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

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 1 ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, solid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 1 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, 1 ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Coarse model

        maxUsedIterations = N;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->solid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 0.01 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 0.01 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        shared_ptr<OutputSpaceMapping> outputSpaceMapping( new OutputSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

        // Create manifold mapping solver
        solver = new SpaceMappingSolver( fineModel, coarseModel, outputSpaceMapping );

        // Monolithic solver
        monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho );

        assert( fluid );
        assert( solid );
        assert( multiLevelFluidSolver );
        assert( multiLevelSolidSolver );
        assert( convergenceMeasures );
        assert( multiLevelFsiSolver );
        assert( postProcessing );
    }

    virtual void TearDown()
    {
        delete solver;
        delete monolithicSolver;
        fineModelFluid.reset();
    }

    MonolithicFsiSolver * monolithicSolver;
    SpaceMappingSolver * solver;
    shared_ptr<TubeFlowFluidSolver> fineModelFluid;
};

TEST_F( OutputSpaceMappingSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( OutputSpaceMappingSolverTest, solveTimeStep )
{
    solver->solveTimeStep();

    ASSERT_LE( solver->fineModel->fsi->nbIter, 10 );
    ASSERT_LE( fineModelFluid->nbRes, 35 );
    ASSERT_LE( fineModelFluid->nbJac, 30 );

    solver->solveTimeStep();

    ASSERT_LE( solver->fineModel->fsi->nbIter, 20 );
    ASSERT_LE( fineModelFluid->nbRes, 100 );
    ASSERT_LE( fineModelFluid->nbJac, 60 );

    solver->solveTimeStep();

    ASSERT_LE( solver->fineModel->fsi->nbIter, 30 );
    ASSERT_LE( fineModelFluid->nbRes, 150 );
    ASSERT_LE( fineModelFluid->nbJac, 90 );

    solver->solveTimeStep();

    ASSERT_LE( solver->fineModel->fsi->nbIter, 40 );
    ASSERT_LE( fineModelFluid->nbRes, 200 );
    ASSERT_LE( fineModelFluid->nbJac, 120 );

    solver->solveTimeStep();

    ASSERT_LE( solver->fineModel->fsi->nbIter, 50 );
    ASSERT_LE( fineModelFluid->nbRes, 250 );
    ASSERT_LE( fineModelFluid->nbJac, 150 );
}
