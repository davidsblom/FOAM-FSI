
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearizedFluidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;

class TubeFlowLinearizedFluidSolverTest : public::testing::Test
{
protected:

  virtual void SetUp()
  {
    double p0 = 0;
    double r0 = 0.2;
    double u0 = 0.1;
    double dt = 0.1;
    double rho = 1.225;
    double L = 1;
    double T = 1;
    double h = 1.0e-3;
    double E = 490;
    int N = 5;

    fluid = new tubeflow::TubeFlowLinearizedFluidSolver( N, p0, r0, u0, rho, E, h, T, dt, L );
  }

  virtual void TearDown()
  {
    delete fluid;
    fluid = NULL;
  }

  TubeFlowLinearizedFluidSolver * fluid;
};

TEST_F( TubeFlowLinearizedFluidSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( TubeFlowLinearizedFluidSolverTest, solve )
{
  fluid->initTimeStep();

  // Initialize variables
  int N = 5;
  fsi::vector p( N ), a( N );
  p.setZero();
  double r0 = 0.2;
  double a0 = M_PI * r0 * r0;
  a.fill( a0 );

  fluid->solve( a, p );
  fluid->finalizeTimeStep();
}
