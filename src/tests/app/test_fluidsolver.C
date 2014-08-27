
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"
#include <math.h>
#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;

class FluidSolverTest : public::testing::Test
{
protected:

  virtual void SetUp()
  {
    double r0 = 0.2;
    double a0 = M_PI * r0 * r0;
    double u0 = 0.1;
    double p0 = 0;
    double dt = 0.1;
    int N = 5;
    double L = 1;
    double T = 1;
    double dx = L / N;
    double rho = 1.225;
    double E = 490;
    double h = 1.0e-3;
    double cmk = std::sqrt( E * h / (2 * rho * r0) );
    double c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
    double kappa = c0 / u0;
    double tau = u0 * dt / L;

    ASSERT_NEAR( tau, 0.01, 1.0e-13 );
    ASSERT_NEAR( kappa, 10, 1.0e-13 );
    ASSERT_TRUE( dx > 0 );

    fluid = new tubeflow::TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho );
  }

  virtual void TearDown()
  {
    delete fluid;
    fluid = NULL;
  }

  TubeFlowFluidSolver * fluid;
};

TEST_F( FluidSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( FluidSolverTest, tau )
{
  ASSERT_NEAR( fluid->tau, 0.01, 1.0e-14 );
}

TEST_F( FluidSolverTest, alpha )
{
  ASSERT_NEAR( fluid->alpha, 0.05983986006, 1.0e-11 );
}

TEST_F( FluidSolverTest, dx )
{
  ASSERT_NEAR( fluid->dx, 0.2, 1.0e-14 );
}

TEST_F( FluidSolverTest, cmk )
{
  ASSERT_NEAR( fluid->cmk, 1, 1.0e-14 );
}

TEST_F( FluidSolverTest, rho )
{
  ASSERT_NEAR( fluid->rho, 1.225, 1.0e-14 );
}

TEST_F( FluidSolverTest, inletVelocityBoundaryCondition )
{
  fluid->init = true;

  double u_in = fluid->evaluateInletVelocityBoundaryCondition();

  ASSERT_NEAR( u_in, 0.1, 1.0e-14 );

  fluid->init = false;
  fluid->initTimeStep();

  u_in = fluid->evaluateInletVelocityBoundaryCondition();

  ASSERT_NEAR( u_in, 0.10000986635, 1.0e-11 );
}

TEST_F( FluidSolverTest, outputPressureBoundaryCondition )
{
  fluid->initTimeStep();

  double p_out = fluid->evaluateOutputPressureBoundaryCondition( 0.5, 1, 1 );

  ASSERT_NEAR( p_out, 0.5, 1.0e-13 );

  p_out = fluid->evaluateOutputPressureBoundaryCondition( 0.5, 2, 2 );

  ASSERT_NEAR( p_out, 0.5, 1.0e-13 );

  p_out = fluid->evaluateOutputPressureBoundaryCondition( 0.5, 2, 1 );

  ASSERT_NEAR( p_out, -0.74599964972, 1.0e-12 );
}

TEST_F( FluidSolverTest, residual )
{
  fluid->initTimeStep();

  int N = 5;
  double r0 = 0.2;
  double a0 = M_PI * r0 * r0;
  double u0 = 0.1;
  double p0 = 0;

  vector R( 2 * N ), x( 2 * N ), a( N ), un( N ), pn( N ), an( N );
  an.fill( a0 );
  pn.fill( p0 );
  un.fill( u0 );
  a = an;
  x.head( N ) = un;
  x.tail( N ) = pn;
  R.setZero();

  // specify un, pn, a, and an by hand

  fluid->evaluateResidual( x, a, un, pn, an, R );

  ASSERT_NEAR( R.norm(), 9.86635786e-06, 1.0e-10 );
}

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

  int inputs() const {
    return m_inputs;
  }

  int values() const {
    return m_values;
  }
};

struct residualFunctor : Functor<double>
{
  residualFunctor(
    TubeFlowFluidSolver * fluid,
    vector * a,
    vector * un,
    vector * pn,
    vector * an
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
  vector * a;
  vector * un;
  vector * pn;
  vector * an;
};

TEST_F( FluidSolverTest, jacobian )
{
  fluid->initTimeStep();

  int N = 5;
  double r0 = 0.2;
  double a0 = M_PI * r0 * r0;
  double u0 = 0.1;
  double p0 = 0;

  vector R( 2 * N ), x( 2 * N ), a( N ), un( N ), pn( N ), an( N );
  matrix J( 2 * N, 2 * N );
  an.fill( a0 );
  pn.fill( p0 );
  un.fill( u0 );
  a = an;
  x.head( N ) = un;
  x.tail( N ) = pn;
  R.setZero();

  fluid->evaluateJacobian( x, a, un, pn, an, J );

  double normJ = J.norm();

  residualFunctor functor( fluid, &a, &un, &pn, &an );
  Eigen::NumericalDiff<residualFunctor, Eigen::Central> numDiff( functor );

  matrix Jnummdiff( 2 * N, 2 * N );
  Jnummdiff.setZero();

  numDiff.df( x, Jnummdiff );

  ASSERT_NEAR( normJ, Jnummdiff.norm(), 1.0e-8 );

  for ( int i = 0; i < J.rows(); i++ )
    for ( int j = 0; j < J.cols(); j++ )
      ASSERT_NEAR( J( i, j ), Jnummdiff( i, j ), 1.0e-7 );
}
