
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "IDWInterpolation.H"
#include "gtest/gtest.h"
#include <math.h>
#include <math.h>

using namespace fsi;

TEST( IDWInterpolationTest, cos )
{
  // Initialize variables

  matrix x( 12, 1 ), y( 12, 1 );

  x << 0, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6;

  for ( int i = 0; i < x.rows(); i++ )
    y( i, 0 ) = std::cos( x( i, 0 ) );

  matrix xnew( 121, 1 ), ynew( 121, 1 );

  for ( int i = 0; i < xnew.rows(); i++ )
    xnew( i, 0 ) = 0.05 * i;

  IDWInterpolation interpolator = IDWInterpolation( 1.8 );

  interpolator.interpolate( x, xnew, y, ynew );

  ASSERT_NEAR( 0.988601249791357, ynew( 1 ), 1.0e-14 );
  ASSERT_NEAR( 0.959479593993784, ynew( 2 ), 1.0e-14 );
  ASSERT_NEAR( 0.915592043333770, ynew( 3 ), 1.0e-14 );
  ASSERT_NEAR( 0.860072479693467, ynew( 4 ), 1.0e-14 );
  ASSERT_NEAR( 0.796738858533311, ynew( 5 ), 1.0e-14 );
  ASSERT_NEAR( 0.729778555289571, ynew( 6 ), 1.0e-14 );
  ASSERT_NEAR( 0.663291723158544, ynew( 7 ), 1.0e-14 );
  ASSERT_NEAR( 0.600889102955681, ynew( 8 ), 1.0e-14 );

  interpolator = IDWInterpolation();

  interpolator.interpolate( x, xnew, y, ynew );

  ASSERT_NEAR( 0.99447731470662892, ynew( 1 ), 1.0e-14 );
}

TEST( IDWInterpolationTest, 3d )
{
  // Check that the RBF function interpolates through the nodes (2D)

  IDWInterpolation idw;

  matrix x, y, ynew;

  x = Eigen::MatrixXd::Random( 50, 3 ).array() * 4 - 2;
  y = Eigen::MatrixXd::Random( 50, 3 ).array() * 4 - 2;

  idw.interpolate( x, x, y, ynew );

  for ( int i = 0; i < y.rows(); i++ )
    for ( int j = 0; j < y.cols(); j++ )
      ASSERT_FALSE( std::isnan( ynew( i, j ) ) );
}
