
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
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"
#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;
using namespace rbf;
using std::shared_ptr;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class ImplicitFsiSolverParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, bool> >
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

    // Parametrized settings
    bool parallel = std::tr1::get<0>( GetParam() );
    int nbReuse = std::tr1::get<1>( GetParam() );
    int extrapolation = std::tr1::get<2>( GetParam() );
    int minIter = std::tr1::get<3>( GetParam() );
    bool convergenceMeasureTraction = std::tr1::get<4>( GetParam() );

    int maxUsedIterations = N;

    if ( parallel )
      maxUsedIterations *= 2;

    ASSERT_NEAR( tau, 0.01, 1.0e-13 );
    ASSERT_NEAR( kappa, 10, 1.0e-13 );
    ASSERT_TRUE( dx > 0 );

    shared_ptr<RBFFunctionInterface> rbfFunction;
    shared_ptr<RBFInterpolation> rbfInterpolator;
    shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
    shared_ptr<RBFCoarsening> rbfInterpToMesh;
    bool coarsening = false;
    int coarseningMinPoints = 200;
    int coarseningMaxPoints = 2000;

    shared_ptr<TubeFlowFluidSolver> fluid( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    shared_ptr<TubeFlowSolidSolver> solid( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

    // Convergence measures
    shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, minIter ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel || convergenceMeasureTraction )
      convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 1, tol ) ) );

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

INSTANTIATE_TEST_CASE_P( testParameters, ImplicitFsiSolverParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Values( 0, 1, 2 ), Values( 1 ), Bool() ) );

TEST_P( ImplicitFsiSolverParametrizedTest, object )
{
  ASSERT_TRUE( true );
}

TEST_P( ImplicitFsiSolverParametrizedTest, initTimeStep )
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

TEST_P( ImplicitFsiSolverParametrizedTest, residual )
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

TEST_P( ImplicitFsiSolverParametrizedTest, iqn_evaluate_residual )
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

TEST_P( ImplicitFsiSolverParametrizedTest, timestep )
{
  solver->solveTimeStep();

  ASSERT_TRUE( solver->fsi->allConverged );
  ASSERT_TRUE( solver->fsi->fluid->isRunning() );
}

TEST_P( ImplicitFsiSolverParametrizedTest, run )
{
  solver->run();

  ASSERT_TRUE( solver->fsi->allConverged );
  ASSERT_FALSE( solver->fsi->fluid->isRunning() );
}

TEST_P( ImplicitFsiSolverParametrizedTest, monolithic )
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

TEST_P( ImplicitFsiSolverParametrizedTest, numberOfColumnsVIQN )
{
  int nbReuse = std::tr1::get<1>( GetParam() );
  int nbIterFirstTimeStep = 0;

  for ( int i = 0; i < 3; i++ )
  {
    solver->initTimeStep();
    solver->solve();

    // Determine the number of columns of the V and W matrices
    int nbResiduals = solver->postProcessing->residuals.size();

    // Include information from previous optimization solves
    for ( std::deque<deque<fsi::vector> >::iterator it = solver->postProcessing->solsList.begin(); it != solver->postProcessing->solsList.end(); ++it )
      nbResiduals += it->size();

    // Include information from previous time steps
    for ( std::deque< std::deque<deque<fsi::vector> > >::iterator solsIterator = solver->postProcessing->solsTimeList.begin(); solsIterator != solver->postProcessing->solsTimeList.end(); ++solsIterator )
      for ( std::deque<deque<fsi::vector> >::iterator it = solsIterator->begin(); it != solsIterator->end(); ++it )
        nbResiduals += it->size();

    if ( i == 0 )
      nbIterFirstTimeStep = solver->fsi->iter;

    if ( nbReuse == 0 )
      ASSERT_EQ( nbResiduals, solver->fsi->iter );

    if ( nbReuse == 1 && i < 2 )
      ASSERT_EQ( nbResiduals, solver->fsi->nbIter );

    if ( nbReuse == 1 && i == 2 )
      ASSERT_EQ( nbResiduals, solver->fsi->nbIter - nbIterFirstTimeStep );

    if ( nbReuse > 1 )
      ASSERT_EQ( nbResiduals, solver->fsi->nbIter );

    solver->finalizeTimeStep();
  }
}

class ImplicitFsiSolverJacobianTest : public::testing::TestWithParam<int>
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
    int N = GetParam();
    double L = 1;
    double T = 10;
    double rho = 1.225;
    double E = 490;
    double h = 1.0e-3;
    double cmk = std::sqrt( E * h / (2 * rho * r0) );

    // Computational settings
    double tol = 1.0e-5;
    int maxIter = 50;
    double initialRelaxation = 1.0e-3;
    double singularityLimit = 1.0e-11;
    bool parallel = false;
    int nbReuse = 0;
    int extrapolation = 0;
    int minIter = 2;
    int reuseInformationStartingFromTimeIndex = 0;
    int maxUsedIterations = N;
    bool scaling = false;
    double beta = 1;
    bool updateJacobian = false;

    shared_ptr<TubeFlowFluidSolver> fluid( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    shared_ptr<TubeFlowSolidSolver> solid( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
    fluid2 = new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho );

    shared_ptr<RBFFunctionInterface> rbfFunction;
    shared_ptr<RBFInterpolation> rbfInterpolator;
    shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
    shared_ptr<RBFCoarsening> rbfInterpToMesh;
    bool coarsening = false;
    int coarseningMinPoints = 200;
    int coarseningMaxPoints = 2000;

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

    // Convergence measures
    shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );
    convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );

    shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );
    shared_ptr<AndersonPostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );
    solver = new ImplicitMultiLevelFsiSolver( fsi, postProcessing );
  }

  virtual void TearDown()
  {
    delete solver;
    delete fluid2;
  }

  ImplicitMultiLevelFsiSolver * solver;
  TubeFlowFluidSolver * fluid2;
};

INSTANTIATE_TEST_CASE_P( jacobian, ImplicitFsiSolverJacobianTest, ::testing::Values( 5, 10, 100, 200 ) );

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() : m_inputs( InputsAtCompileTime ), m_values( ValuesAtCompileTime ) {}

  Functor(
    int inputs,
    int values
    ) : m_inputs( inputs ), m_values( values ) {}

  int inputs() const
  {
    return m_inputs;
  }

  int values() const
  {
    return m_values;
  }
};

struct residualFunctor : Functor<double>
{
  residualFunctor(
    TubeFlowFluidSolver * fluid,
    fsi::vector * a,
    fsi::vector * un,
    fsi::vector * pn,
    fsi::vector * an
    )
    :
    Functor<double>( 2 * fluid->N, 2 * fluid->N ),
    fluid( fluid ),
    a( a ),
    un( un ),
    pn( pn ),
    an( an )
  {}

  int operator()(
    Eigen::VectorXd & x,
    Eigen::VectorXd & fvec
    ) const
  {
    fluid->evaluateResidual( x, *a, *un, *pn, *an, fvec );
    return 0;
  }

  TubeFlowFluidSolver * fluid;
  fsi::vector * a;
  fsi::vector * un;
  fsi::vector * pn;
  fsi::vector * an;
};


TEST_P( ImplicitFsiSolverJacobianTest, jacobian )
{
  // Verify the Jacobian matrix of the fluid solver at each time step

  for ( int i = 0; i < 100; i++ )
  {
    int N = GetParam();

    solver->solveTimeStep();

    ASSERT_TRUE( solver->fsi->allConverged );

    if ( i < 99 )
      ASSERT_TRUE( solver->fsi->fluid->isRunning() );
    else
      ASSERT_FALSE( solver->fsi->fluid->isRunning() );

    fluid2->initTimeStep();

    matrix p( N, 1 );

    fluid2->solve( solver->fsi->solid->data, p );

    matrix J( 2 * N, 2 * N );
    fsi::vector R( 2 * N );
    fsi::vector x( 2 * N );
    x.head( N ) = fluid2->u;
    x.tail( N ) = fluid2->p;

    fluid2->evaluateJacobian( x, fluid2->a, fluid2->un, fluid2->pn, fluid2->an, J );

    double normJ = J.norm();

    residualFunctor functor( fluid2, &(fluid2->a), &(fluid2->un), &(fluid2->pn), &(fluid2->an) );

    // test the functor
    functor( x, R );

    Eigen::NumericalDiff<residualFunctor, Eigen::Central> numDiff( functor );

    // Eigen::NumericalDiff<residualFunctor, Eigen::FourthOrder> numDiff(functor);

    matrix Jnummdiff( 2 * N, 2 * N );
    Jnummdiff.setZero();

    numDiff.df( x, Jnummdiff );

    ASSERT_NEAR( normJ, Jnummdiff.norm(), 1.0e-6 );

    for ( int i = 0; i < J.rows(); i++ )
      for ( int j = 0; j < J.cols(); j++ )
        ASSERT_NEAR( J( i, j ), Jnummdiff( i, j ), 1.0e-6 );

    fluid2->finalizeTimeStep();
  }
}
