
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "BroydenPostProcessing.H"
#include "ManifoldMapping.H"
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

class ManifoldMappingSolverParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, int, bool, bool> >
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
    bool broyden = std::tr1::get<5>( GetParam() );
    bool convergenceMeasureTraction = std::tr1::get<6>( GetParam() );

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
    shared_ptr<PostProcessing> postProcessing;
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

    if ( parallel || convergenceMeasureTraction )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    if ( broyden )
      postProcessing = shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    if ( !broyden )
      postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

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

    if ( parallel || convergenceMeasureTraction )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

    multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    if ( broyden )
      postProcessing = shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    if ( !broyden )
      postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    // Create manifold mapping object

    shared_ptr<ManifoldMapping> manifoldMapping( new ManifoldMapping( fineModel, coarseModel, maxIter, singularityLimit, nbReuse, reuseInformationStartingFromTimeIndex ) );

    // Create manifold mapping solver
    solver = new SpaceMappingSolver( fineModel, coarseModel, manifoldMapping );

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

INSTANTIATE_TEST_CASE_P( testParameters, ManifoldMappingSolverParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Values( 0 ), Values( 3 ), Values( 10, 20 ), Bool(), Bool() ) );

TEST_P( ManifoldMappingSolverParametrizedTest, object )
{
  ASSERT_TRUE( true );
}

TEST_P( ManifoldMappingSolverParametrizedTest, run )
{
  solver->run();

  ASSERT_TRUE( solver->fineModel->fsi->allConverged );
  ASSERT_FALSE( solver->fineModel->fsi->fluid->isRunning() );
  ASSERT_GT( solver->coarseModel->fsi->nbIter, solver->fineModel->fsi->nbIter );

  int couplingGridSize = std::tr1::get<4>( GetParam() );
  bool convergenceMeasureTraction = std::tr1::get<6>( GetParam() );

  if ( couplingGridSize == 10 && !convergenceMeasureTraction )
    ASSERT_EQ( solver->fineModel->fsi->nbIter, 100 );

  ASSERT_GT( solver->timeElapsed(), 0 );
  ASSERT_EQ( solver->spaceMapping->timeIndex, 100 );
  ASSERT_EQ( solver->fineModel->postProcessing->timeIndex, 100 );
  ASSERT_EQ( solver->coarseModel->postProcessing->timeIndex, 100 );
}

TEST_P( ManifoldMappingSolverParametrizedTest, monolithic )
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
    ASSERT_TRUE( solver->coarseModel->fsi->allConverged );
    ASSERT_NEAR( solver->fineModel->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
    ASSERT_NEAR( solver->fineModel->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
    ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
    ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );
  }
}

TEST_P( ManifoldMappingSolverParametrizedTest, postProcessingCoarseModel )
{
  int nbReuse = std::tr1::get<1>( GetParam() );
  int nbIterFirstTimeStep = 0;
  int nbIterPreviousTimestep = 0;

  for ( int i = 0; i < 3; i++ )
  {
    solver->initTimeStep();
    solver->solve();

    // Determine the number of columns of the V and W matrices
    int nbResiduals = solver->coarseModel->postProcessing->residuals.size();

    // Include information from previous optimization solves
    for ( std::deque<deque<fsi::vector> >::iterator it = solver->coarseModel->postProcessing->solsList.begin(); it != solver->coarseModel->postProcessing->solsList.end(); ++it )
      nbResiduals += it->size();

    // Include information from previous time steps
    for ( std::deque< std::deque<deque<fsi::vector> > >::iterator solsIterator = solver->coarseModel->postProcessing->solsTimeList.begin(); solsIterator != solver->coarseModel->postProcessing->solsTimeList.end(); ++solsIterator )
      for ( std::deque<deque<fsi::vector> >::iterator it = solsIterator->begin(); it != solsIterator->end(); ++it )
        nbResiduals += it->size();

    if ( i == 0 )
      nbIterFirstTimeStep = solver->coarseModel->fsi->iter;

    if ( nbReuse == 0 )
      ASSERT_EQ( nbResiduals, solver->coarseModel->fsi->nbIter - solver->fineModel->fsi->iter - nbIterPreviousTimestep );

    if ( nbReuse == 1 && i < 2 )
      ASSERT_EQ( nbResiduals, solver->coarseModel->fsi->nbIter - solver->fineModel->fsi->nbIter );

    if ( nbReuse == 1 && i == 2 )
      ASSERT_EQ( nbResiduals, solver->coarseModel->fsi->nbIter - solver->fineModel->fsi->nbIter - nbIterFirstTimeStep );

    if ( nbReuse > 1 )
      ASSERT_EQ( nbResiduals, solver->coarseModel->fsi->nbIter - solver->fineModel->fsi->nbIter );

    nbIterPreviousTimestep = solver->coarseModel->fsi->nbIter;

    if ( i == 0 )
      nbIterFirstTimeStep = solver->coarseModel->fsi->nbIter - solver->fineModel->fsi->nbIter;

    solver->finalizeTimeStep();
  }
}

TEST_P( ManifoldMappingSolverParametrizedTest, postProcessingFineModel )
{
  int nbReuse = std::tr1::get<1>( GetParam() );

  for ( int i = 0; i < 3; i++ )
  {
    solver->initTimeStep();
    solver->solve();

    // Determine the number of columns of the V and W matrices
    int nbResiduals = solver->spaceMapping->fineResiduals.cols();

    // Include information from previous time steps
    for ( std::deque<matrix>::iterator it = solver->spaceMapping->fineResidualsList.begin(); it != solver->spaceMapping->fineResidualsList.end(); ++it )
      nbResiduals += it->cols();

    if ( nbReuse == 0 )
      ASSERT_EQ( nbResiduals, 0 );

    if ( solver->fineModel->fsi->iter > 1 )
    {
      if ( nbReuse == 1 )
        ASSERT_EQ( nbResiduals, solver->fineModel->fsi->iter );

      if ( nbReuse > 1 )
        ASSERT_EQ( nbResiduals, solver->fineModel->fsi->nbIter );
    }

    solver->finalizeTimeStep();
  }
}


class ManifoldMappingSolverTest : public::testing::Test
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

    // Parametrized settings
    bool parallel = false;
    int nbReuse = 0;
    int extrapolation = 0;
    int minIter = 2;
    int couplingGridSize = 100;
    bool scaling = false;

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
    shared_ptr<PostProcessing> postProcessing;
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
    postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

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

    postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    // Create manifold mapping object

    shared_ptr<ManifoldMapping> manifoldMapping( new ManifoldMapping( fineModel, coarseModel, maxIter, singularityLimit, nbReuse, reuseInformationStartingFromTimeIndex ) );

    // Create manifold mapping solver
    solver = new SpaceMappingSolver( fineModel, coarseModel, manifoldMapping );

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

TEST_F( ManifoldMappingSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( ManifoldMappingSolverTest, solveTimeStep )
{
  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 9 );
  ASSERT_EQ( fineModelFluid->nbRes, 45 );
  ASSERT_EQ( fineModelFluid->nbJac, 27 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 18 );
  ASSERT_EQ( fineModelFluid->nbRes, 90 );
  ASSERT_EQ( fineModelFluid->nbJac, 54 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 26 );
  ASSERT_EQ( fineModelFluid->nbRes, 130 );
  ASSERT_EQ( fineModelFluid->nbJac, 78 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 34 );
  ASSERT_EQ( fineModelFluid->nbRes, 170 );
  ASSERT_EQ( fineModelFluid->nbJac, 102 );

  solver->solveTimeStep();

  ASSERT_EQ( solver->fineModel->fsi->nbIter, 41 );
  ASSERT_EQ( fineModelFluid->nbRes, 205 );
  ASSERT_EQ( fineModelFluid->nbJac, 123 );
}
