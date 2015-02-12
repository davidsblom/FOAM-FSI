
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MonolithicFsiSolver.H"
#include "gtest/gtest.h"
using namespace tubeflow;

class MonolithicFsiSolverTest : public::testing::Test
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
    double rho = 1.225;
    double E = 490;
    double h = 1.0e-3;
    double cmk = std::sqrt( E * h / (2 * rho * r0) );

    fsi = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, N, L, T, rho );
  }

  virtual void TearDown()
  {
    delete fsi;
  }

  MonolithicFsiSolver * fsi;
};

TEST_F( MonolithicFsiSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( MonolithicFsiSolverTest, run )
{
  fsi->run();

  double tol = 1.0e-5;

  fsi::vector pn( 5 ), un( 5 ), an( 5 );
  pn << -1.04294792e-06, 1.11864827e-04, 2.24772603e-04, 4.52098904e-04, 7.13923087e-04;
  un << 0.1, 0.09999385, 0.10015137, 0.1003671, 0.10058284;
  an << 0.1256636, 0.12567518, 0.12568677, 0.1257101, 0.12573697;

  for ( int i = 0; i < 5; i++ )
  {
    ASSERT_NEAR( fsi->pn( i ), pn( i ), tol );
    ASSERT_NEAR( fsi->un( i ), un( i ), tol );
    ASSERT_NEAR( fsi->an( i ), an( i ), tol );
  }
}
