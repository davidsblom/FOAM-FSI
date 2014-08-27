
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"
#include "gtest/gtest.h"

using namespace fsi;

TEST( RBFInterpolationTest, oneD_small )
{
  // Initialize variables

  int n_A = 7;
  int n_B = 6;

  matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

  // Location of original values
  x << 0, 1, 2, 3, 4, 5, 6;
  xnew << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;

  // Values
  y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

  RBFInterpolation rbf( x, xnew, "thin_plate", 1.0 );

  rbf.interpolate( y, ynew );

  ASSERT_NEAR( 0.44626242465329269, ynew( 0 ), 1.0e-14 );
  ASSERT_NEAR( 0.98177874060900239, ynew( 1 ), 1.0e-14 );
  ASSERT_NEAR( 0.58634017875503497, ynew( 2 ), 1.0e-14 );
  ASSERT_NEAR( -0.34305440473301307, ynew( 3 ), 1.0e-14 );
  ASSERT_NEAR( -0.96344050977208218, ynew( 4 ), 1.0e-14 );
  ASSERT_NEAR( -0.65619975149611487, ynew( 5 ), 1.0e-14 );
  ASSERT_EQ( 1, ynew.cols() );
  ASSERT_EQ( 6, ynew.rows() );
  ASSERT_EQ( xnew.cols(), ynew.cols() );
  ASSERT_EQ( xnew.rows(), ynew.rows() );

  // Verify consistent interpolation. Rowsum of Hhat is one for a consistent
  // interpolation.
  vector rowsum = rbf.Hhat.rowwise().sum();

  for ( int i = 0; i < rowsum.rows(); i++ )
    ASSERT_NEAR( rowsum( i ), 1, 1.0e-14 );
}

TEST( RBFInterpolationTest, oneD_smallSamePosition )
{
  // Initialize variables

  int n_A = 7;
  int n_B = n_A;

  matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

  x << 0, 1, 2, 3, 4, 5, 6;
  xnew = x;

  y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

  RBFInterpolation rbf( x, xnew, "thin_plate", 1.0 );

  rbf.interpolate( y, ynew );

  for ( int i = 0; i < ynew.rows(); i++ )
    ASSERT_NEAR( y( i ), ynew( i ), 1.0e-14 );

  ASSERT_EQ( 1, ynew.cols() );
  ASSERT_EQ( 7, ynew.rows() );
  ASSERT_EQ( xnew.cols(), ynew.cols() );
  ASSERT_EQ( xnew.rows(), ynew.rows() );

  // Verify consistent interpolation. Rowsum of Hhat is one for a consistent
  // interpolation.
  vector rowsum = rbf.Hhat.rowwise().sum();

  for ( int i = 0; i < rowsum.rows(); i++ )
    ASSERT_NEAR( rowsum( i ), 1, 1.0e-14 );
}

TEST( RBFInterpolationTest, rbf1d )
{
  // Check that the RBF function interpolates through the nodes (1D)

  RBFInterpolation rbf;

  matrix x( 9, 1 ), y( 9, 1 ), ynew;
  x.col( 0 ).setLinSpaced( 9, 0, 10 );

  y.col( 0 ) = x.array().sin();

  rbf.compute( x, x );
  rbf.interpolate( y, ynew );

  for ( int i = 0; i < y.rows(); i++ )
    for ( int j = 0; j < y.cols(); j++ )
      ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, rbf2d )
{
  // Check that the RBF function interpolates through the nodes (2D)

  RBFInterpolation rbf;

  matrix x, y, ynew;

  x = Eigen::MatrixXd::Random( 50, 2 ).array() * 4 - 2;
  y = Eigen::MatrixXd::Random( 50, 2 ).array() * 4 - 2;

  rbf.compute( x, x );
  rbf.interpolate( y, ynew );

  for ( int i = 0; i < y.rows(); i++ )
    for ( int j = 0; j < y.cols(); j++ )
      ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, rbf3d )
{
  // Check that the RBF function interpolates through the nodes (2D)

  RBFInterpolation rbf;

  matrix x, y, ynew;

  x = Eigen::MatrixXd::Random( 50, 3 ).array() * 4 - 2;
  y = Eigen::MatrixXd::Random( 50, 3 ).array() * 4 - 2;

  rbf.compute( x, x );
  rbf.interpolate( y, ynew );

  for ( int i = 0; i < y.rows(); i++ )
    for ( int j = 0; j < y.cols(); j++ )
      ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, rbf1d_regularity )
{
  // Check that the RBF function approximates a smooth function well
  // away from the nodes

  matrix x( 9, 1 ), y( 9, 1 ), xnew( 100, 1 ), ynew;
  x.col( 0 ).setLinSpaced( 9, 0, 10 );
  y.col( 0 ) = x.array().sin();
  xnew.col( 0 ).setLinSpaced( 100, 0, 10 );

  RBFInterpolation rbf;

  rbf.compute( x, xnew );
  rbf.interpolate( y, ynew );

  for ( int i = 0; i < ynew.rows(); i++ )
    ASSERT_NEAR( ynew( i, 0 ), std::sin( xnew( i, 0 ) ), 1.0e-1 );
}
