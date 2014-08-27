
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
#include "BroydenPostProcessing.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class OutputSpaceMappingSolverParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, int> >
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
    int reuseInformationStartingFromTimeIndex = 0;

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
    shared_ptr<BroydenPostProcessing> postProcessing;
    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

    int maxUsedIterations = couplingGridSize;

    if ( parallel )
      maxUsedIterations *= 2;

    // Fine model

    fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );
    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, 0, 1 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, solid, 1, 1 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, 1 ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
    postProcessing = shared_ptr<BroydenPostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    // Coarse model

    maxUsedIterations = N;

    if ( parallel )
      maxUsedIterations *= 2;

    fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, 0, 0 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->solid, 1, 0 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    postProcessing = shared_ptr<BroydenPostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );
    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    // Create manifold mapping object

    shared_ptr<OutputSpaceMapping> outputSpaceMapping( new OutputSpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

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

INSTANTIATE_TEST_CASE_P( testParameters, OutputSpaceMappingSolverParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Values( 0 ), Values( 3 ), Values( 10, 20 ) ) );

TEST_P( OutputSpaceMappingSolverParametrizedTest, object )
{
  ASSERT_TRUE( true );
}

TEST_P( OutputSpaceMappingSolverParametrizedTest, run )
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

class OutputSpaceMappingSolverTest : public::testing::Test
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
    int reuseInformationStartingFromTimeIndex = 0;

    // Parametrized settings
    bool parallel = false;
    int nbReuse = 0;
    int extrapolation = 0;
    int minIter = 5;
    int couplingGridSize = 100;

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
    shared_ptr<BroydenPostProcessing> postProcessing;
    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

    // Fine model

    int maxUsedIterations = couplingGridSize;

    fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );
    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, 1 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, solid, 1, 1 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, 1 ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
    postProcessing = shared_ptr<BroydenPostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    // Coarse model

    maxUsedIterations = N;

    fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, 0, 0 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->solid, 1, 0 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    postProcessing = shared_ptr<BroydenPostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    // Create manifold mapping object

    shared_ptr<OutputSpaceMapping> outputSpaceMapping( new OutputSpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

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
  shared_ptr<TubeFlowFluidSolver> fineModelFluid;
};

TEST_F( OutputSpaceMappingSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( OutputSpaceMappingSolverTest, solveTimeStep )
{
  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 11 );
  ASSERT_EQ( fineModelFluid->nbRes, 55 );
  ASSERT_EQ( fineModelFluid->nbJac, 33 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 23 );
  ASSERT_EQ( fineModelFluid->nbRes, 115 );
  ASSERT_EQ( fineModelFluid->nbJac, 69 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 33 );
  ASSERT_EQ( fineModelFluid->nbRes, 165 );
  ASSERT_EQ( fineModelFluid->nbJac, 99 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 43 );
  ASSERT_EQ( fineModelFluid->nbRes, 215 );
  ASSERT_EQ( fineModelFluid->nbJac, 129 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 54 );
  ASSERT_EQ( fineModelFluid->nbRes, 270 );
  ASSERT_EQ( fineModelFluid->nbJac, 162 );
}
