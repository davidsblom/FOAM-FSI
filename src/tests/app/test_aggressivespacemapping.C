
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AndersonPostProcessing.H"
#include "AggressiveSpaceMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
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

class AggressiveSpaceMappingSolverParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, int> >
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
        double tol = 1.0e-5;
        int maxIter = 500;
        double initialRelaxation = 1.0e-3;
        double singularityLimit = 1.0e-11;
        int reuseInformationStartingFromTimeIndex = 0;
        bool scaling = false;
        double beta = 1;
        bool updateJacobian = false;

        // Parametrized settings
        bool parallel = std::tr1::get<0>( GetParam() );
        int nbReuse = std::tr1::get<1>( GetParam() );
        int extrapolation = std::tr1::get<2>( GetParam() );
        int minIter = std::tr1::get<3>( GetParam() );
        int couplingGridSize = std::tr1::get<4>( GetParam() );

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

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, 1 ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

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

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, 1.0e-2 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, 1.0e-2 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create aggressive space mapping object

        maxUsedIterations = couplingGridSize;

        if ( parallel )
            maxUsedIterations *= 2;

        shared_ptr<AggressiveSpaceMapping> aggressiveSpaceMapping( new AggressiveSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

        // Create manifold mapping solver
        solver = new SpaceMappingSolver( fineModel, coarseModel, aggressiveSpaceMapping );

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

INSTANTIATE_TEST_CASE_P( testParameters, AggressiveSpaceMappingSolverParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Values( 0 ), Values( 3 ), Values( 10, 20 ) ) );

TEST_P( AggressiveSpaceMappingSolverParametrizedTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( AggressiveSpaceMappingSolverParametrizedTest, run )
{
    solver->run();
    monolithicSolver->run();

    double tol = 1.0e-5;
    ASSERT_TRUE( solver->fineModel->fsi->allConverged );
    ASSERT_TRUE( solver->coarseModel->fsi->allConverged );
    ASSERT_NEAR( solver->fineModel->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
    ASSERT_NEAR( solver->fineModel->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
    ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
    ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );
}

TEST_P( AggressiveSpaceMappingSolverParametrizedTest, reuse )
{
    solver->run();

    bool parallel = std::tr1::get<0>( GetParam() );
    int nbReuse = std::tr1::get<1>( GetParam() );
    int extrapolation = std::tr1::get<2>( GetParam() );
    int minIter = std::tr1::get<3>( GetParam() );
    int couplingGridSize = std::tr1::get<4>( GetParam() );

    if ( !parallel && couplingGridSize == 20 && nbReuse == 0 && extrapolation == 0 && minIter == 3 )
    {
        ASSERT_EQ( solver->fineModel->fsi->nbIter, 729 );
        ASSERT_EQ( solver->coarseModel->fsi->nbIter, 4691 );
    }

    if ( !parallel && couplingGridSize == 20 && nbReuse == 1 && extrapolation == 0 && minIter == 3 )
    {
        ASSERT_EQ( solver->fineModel->fsi->nbIter, 624 );
        ASSERT_EQ( solver->coarseModel->fsi->nbIter, 3612 );
    }

    if ( !parallel && couplingGridSize == 20 && nbReuse == 4 && extrapolation == 0 && minIter == 3 )
    {
        ASSERT_EQ( solver->fineModel->fsi->nbIter, 543 );
        ASSERT_EQ( solver->coarseModel->fsi->nbIter, 3101 );
    }
}

class AggressiveSpaceMappingSolverTest : public::testing::Test
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
        double tol = 1.0e-5;
        int maxIter = 500;
        double initialRelaxation = 1.0e-3;
        double singularityLimit = 1.0e-11;
        int reuseInformationStartingFromTimeIndex = 0;
        bool updateJacobian = false;

        // Parametrized settings
        bool parallel = false;
        int nbReuse = 0;
        int extrapolation = 0;
        int minIter = 2;
        int couplingGridSize = 100;
        bool scaling = false;
        double beta = 1;

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

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, 1 ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

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

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        shared_ptr<AggressiveSpaceMapping> aggressiveSpaceMapping( new AggressiveSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

        // Create manifold mapping solver
        solver = new SpaceMappingSolver( fineModel, coarseModel, aggressiveSpaceMapping );

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
    shared_ptr<TubeFlowFluidSolver> fineModelFluid;
};

TEST_F( AggressiveSpaceMappingSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( AggressiveSpaceMappingSolverTest, solveTimeStep )
{
    solver->solveTimeStep();

    ASSERT_EQ( solver->fineModel->fsi->nbIter, 9 );
    ASSERT_EQ( fineModelFluid->nbRes, 45 );
    ASSERT_EQ( fineModelFluid->nbJac, 27 );

    solver->solveTimeStep();

    ASSERT_EQ( solver->fineModel->fsi->nbIter, 17 );
    ASSERT_EQ( fineModelFluid->nbRes, 85 );
    ASSERT_EQ( fineModelFluid->nbJac, 51 );

    solver->solveTimeStep();

    ASSERT_EQ( solver->fineModel->fsi->nbIter, 26 );
    ASSERT_EQ( fineModelFluid->nbRes, 130 );
    ASSERT_EQ( fineModelFluid->nbJac, 78 );

    solver->solveTimeStep();

    ASSERT_EQ( solver->fineModel->fsi->nbIter, 35 );
    ASSERT_EQ( fineModelFluid->nbRes, 175 );
    ASSERT_EQ( fineModelFluid->nbJac, 105 );

    solver->solveTimeStep();

    ASSERT_EQ( solver->fineModel->fsi->nbIter, 44 );
    ASSERT_EQ( solver->coarseModel->fsi->nbIter, 232 );
    ASSERT_EQ( fineModelFluid->nbRes, 220 );
    ASSERT_EQ( fineModelFluid->nbJac, 132 );
}
