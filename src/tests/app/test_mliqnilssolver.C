
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MLIQNILSSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class MLIQNILSSolverParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, int> >
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
    shared_ptr<TubeFlowFluidSolver> fineModelFluid;
    shared_ptr<TubeFlowSolidSolver> fineModelSolid;
    shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
    shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
    shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
    shared_ptr<IQNILSPostProcessing> postProcessing;
    shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

    models = shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >  > ( new std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );

    fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
    fineModelSolid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

    // Coarse model

    int maxUsedIterations = N;

    if ( parallel )
      maxUsedIterations *= 2;

    fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, 0, 0 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, 1, 0 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    postProcessing = shared_ptr<IQNILSPostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    assert( fluid );
    assert( solid );
    assert( multiLevelFluidSolver );
    assert( multiLevelSolidSolver );
    assert( convergenceMeasures );
    assert( multiLevelFsiSolver );
    assert( postProcessing );

    fluid.reset();
    solid.reset();
    multiLevelFluidSolver.reset();
    multiLevelSolidSolver.reset();
    convergenceMeasures.reset();
    multiLevelFsiSolver.reset();
    postProcessing.reset();

    // Fine model

    maxUsedIterations = couplingGridSize;

    if ( parallel )
      maxUsedIterations *= 2;

    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, 1 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, 1, 1 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, 1 ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
    postProcessing = shared_ptr<IQNILSPostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    models->push_back( coarseModel );
    models->push_back( fineModel );

    // Create manifold mapping solver
    solver = new MLIQNILSSolver( models );

    // Monolithic solver
    monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho );

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
  MLIQNILSSolver * solver;
};

INSTANTIATE_TEST_CASE_P( testParameters, MLIQNILSSolverParametrizedTest, ::testing::Combine( Bool(), Values( 0, 2 ), Values( 0 ), Values( 3 ), Values( 10, 20 ) ) );

TEST_P( MLIQNILSSolverParametrizedTest, object )
{
  ASSERT_TRUE( true );
}

TEST_P( MLIQNILSSolverParametrizedTest, run )
{
  solver->run();

  ASSERT_TRUE( solver->fineModel->fsi->allConverged );
  ASSERT_FALSE( solver->fineModel->fsi->fluid->isRunning() );

  int couplingGridSize = std::tr1::get<4>( GetParam() );

  if ( couplingGridSize == 10 )
    ASSERT_EQ( solver->fineModel->fsi->nbIter, 100 );
}

TEST_P( MLIQNILSSolverParametrizedTest, monolithic )
{
  for ( int i = 0; i < 100; i++ )
  {
    solver->solveTimeStep();
    monolithicSolver->solveTimeStep();

    double tol = 1.0e-5;

    if ( i < 99 )
      ASSERT_TRUE( solver->fineModel->fsi->fluid->isRunning() );
    else
      ASSERT_FALSE( solver->fineModel->fsi->fluid->isRunning() );

    ASSERT_TRUE( solver->fineModel->fsi->allConverged );
    ASSERT_TRUE( solver->models->at( 0 )->fsi->allConverged );
    ASSERT_NEAR( solver->fineModel->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
    ASSERT_NEAR( solver->fineModel->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
    ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
    ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );

    // Verify that the coarse models are synchronized with the fine model
    for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = solver->models->begin(); it != solver->models->end(); ++it )
    {
      shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;
      ASSERT_NEAR( model->fsi->x.norm(), solver->fineModel->fsi->x.norm(), 1.0e-13 );
    }
  }
}

TEST_P( MLIQNILSSolverParametrizedTest, postProcessingTimeStep )
{
  for ( int i = 0; i < 1; i++ )
  {
    solver->initTimeStep();
    solver->solve();

    // Verify the number of columns of the V and W matrices

    int nbResiduals = solver->fineModel->postProcessing->residuals.size();

    // Include information from previous optimization solves
    for ( std::deque<deque<fsi::vector> >::iterator it = solver->fineModel->postProcessing->solsList.begin(); it != solver->fineModel->postProcessing->solsList.end(); ++it )
      nbResiduals += it->size();

    int nbIter = 0;

    for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = solver->models->begin(); it != solver->models->end(); ++it )
    {
      shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;

      if ( model->fsi->nbIter > 1 )
        nbIter += model->fsi->nbIter;
    }

    ASSERT_EQ( nbIter, nbResiduals );
  }
}

TEST_P( MLIQNILSSolverParametrizedTest, postProcessingMultipleTimeSteps )
{
  int nbReuse = std::tr1::get<1>( GetParam() );
  int nbIter = 0;
  int nbIterFirstTimeStep = 0;

  for ( int i = 0; i < 4; i++ )
  {
    solver->initTimeStep();
    solver->solve();

    // Verify the number of columns of the V and W matrices

    int nbResiduals = solver->fineModel->postProcessing->residuals.size();

    // Include information from previous optimization solves
    for ( std::deque<deque<fsi::vector> >::iterator it = solver->fineModel->postProcessing->solsList.begin(); it != solver->fineModel->postProcessing->solsList.end(); ++it )
      nbResiduals += it->size();

    // Include information from previous time steps
    for ( std::deque< std::deque<deque<fsi::vector> > >::iterator solsIterator = solver->fineModel->postProcessing->solsTimeList.begin(); solsIterator != solver->fineModel->postProcessing->solsTimeList.end(); ++solsIterator )
      for ( std::deque<deque<fsi::vector> >::iterator it = solsIterator->begin(); it != solsIterator->end(); ++it )
        nbResiduals += it->size();

    if ( nbReuse == 0 )
      nbIter = 0;

    for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = solver->models->begin(); it != solver->models->end(); ++it )
    {
      shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;

      if ( model->fsi->iter > 1 )
        nbIter += model->fsi->iter;
    }

    if ( nbReuse == 2 && i == 3 )
      nbIter -= nbIterFirstTimeStep;

    if ( i == 0 )
      nbIterFirstTimeStep = nbIter;

    ASSERT_EQ( nbIter, nbResiduals );

    solver->finalizeTimeStep();
  }
}

class MLIQNILSSolverTest : public::testing::Test
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

    // Parametrized settings
    bool parallel = false;
    int nbReuse = 0;
    int extrapolation = 0;
    int minIter = 1;
    int couplingGridSize = 100;
    int reuseInformationStartingFromTimeIndex = 0;
    bool scaling = false;

    ASSERT_NEAR( tau, 0.01, 1.0e-13 );
    ASSERT_NEAR( kappa, 10, 1.0e-13 );
    ASSERT_TRUE( dx > 0 );

    // Initialize variables
    shared_ptr<TubeFlowSolidSolver> solid;
    shared_ptr<TubeFlowFluidSolver> fluid;
    shared_ptr<TubeFlowSolidSolver> fineModelSolid;
    shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
    shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
    shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
    shared_ptr<IQNILSPostProcessing> postProcessing;
    shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

    models = shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >  > ( new std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );

    fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
    fineModelSolid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

    // Coarse model

    int maxUsedIterations = N;

    if ( parallel )
      maxUsedIterations *= 2;

    fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, 0, 0 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, 1, 0 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    postProcessing = shared_ptr<IQNILSPostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    assert( fluid );
    assert( solid );
    assert( multiLevelFluidSolver );
    assert( multiLevelSolidSolver );
    assert( convergenceMeasures );
    assert( multiLevelFsiSolver );
    assert( postProcessing );

    fluid.reset();
    solid.reset();
    multiLevelFluidSolver.reset();
    multiLevelSolidSolver.reset();
    convergenceMeasures.reset();
    multiLevelFsiSolver.reset();
    postProcessing.reset();

    // Fine model

    maxUsedIterations = couplingGridSize;

    if ( parallel )
      maxUsedIterations *= 2;

    multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, 1 ) );
    multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, 1, 1 ) );

    // Convergence measures
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, 1 ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
    postProcessing = shared_ptr<IQNILSPostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    models->push_back( coarseModel );
    models->push_back( fineModel );

    // Create manifold mapping solver
    solver = new MLIQNILSSolver( models );

    // Monolithic solver
    monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho );

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
  MLIQNILSSolver * solver;
  shared_ptr<TubeFlowFluidSolver> fineModelFluid;
};

TEST_F( MLIQNILSSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( MLIQNILSSolverTest, solveTimeStep )
{
  solver->solveTimeStep();

  ASSERT_EQ( solver->models->back()->fsi->nbIter, 8 );
  ASSERT_EQ( fineModelFluid->nbRes, 40 );
  ASSERT_EQ( fineModelFluid->nbJac, 24 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->models->back()->fsi->nbIter, 18 );
  ASSERT_EQ( fineModelFluid->nbRes, 90 );
  ASSERT_EQ( fineModelFluid->nbJac, 54 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->models->back()->fsi->nbIter, 28 );
  ASSERT_EQ( fineModelFluid->nbRes, 140 );
  ASSERT_EQ( fineModelFluid->nbJac, 84 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->models->back()->fsi->nbIter, 37 );
  ASSERT_EQ( fineModelFluid->nbRes, 185 );
  ASSERT_EQ( fineModelFluid->nbJac, 111 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->models->back()->fsi->nbIter, 46 );
  ASSERT_EQ( fineModelFluid->nbRes, 230 );
  ASSERT_EQ( fineModelFluid->nbJac, 138 );
}
