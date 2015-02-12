
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;

class SolidSolverTest : public::testing::Test
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
    ASSERT_NEAR( cmk, 1, 1.0e-15 );
    ASSERT_TRUE( dx > 0 );

    solid = new tubeflow::TubeFlowSolidSolver( a0, cmk, p0, rho, L, N );
  }

  virtual void TearDown()
  {
    delete solid;
    solid = NULL;
  }

  TubeFlowSolidSolver * solid;
};

TEST_F( SolidSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( SolidSolverTest, timestep )
{
  solid->initTimeStep();

  // Initialize variables
  double r0 = 0.2;
  double a0 = M_PI * r0 * r0;
  int N = 5;
  fsi::vector p( N ), a( N );
  p.setZero();
  a.setZero();

  solid->solve( p, a );

  for ( int i = 0; i < a.rows(); i++ )
    ASSERT_NEAR( a( i ), a0, 1.0e-15 );

  p.fill( 1 );

  solid->solve( p, a );

  for ( int i = 0; i < a.rows(); i++ )
    ASSERT_NEAR( a( i ), 0.35876166284, 1.0e-11 );

  solid->finalizeTimeStep();
}

TEST_F( SolidSolverTest, po )
{
  solid->p0 = 2;
  solid->initTimeStep();

  // Initialize variables
  int N = 5;
  fsi::vector p( N ), a( N );
  p.setZero();
  a.setZero();

  solid->solve( p, a );

  for ( int i = 0; i < a.rows(); i++ )
    ASSERT_NEAR( a( i ), 0.00423938367, 1.0e-11 );

  solid->finalizeTimeStep();
}
