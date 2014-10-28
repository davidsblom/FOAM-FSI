
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "IDWInterpolation.H"
#include <math.h>

namespace fsi
{
  IDWInterpolation::IDWInterpolation() :
    c( 2 )
  {}

  IDWInterpolation::IDWInterpolation( double c ) :
    c( c )
  {}

  IDWInterpolation::~IDWInterpolation()
  {}

  void IDWInterpolation::interpolate(
    const matrix & positions,
    const matrix & positionsInterpolation,
    const matrix & values,
    matrix & valuesInterpolation
    )
  {
    // Initialize variables

    int n = positions.rows();
    int m = positionsInterpolation.rows();

    valuesInterpolation.resize( positionsInterpolation.rows(), values.cols() );
    valuesInterpolation.setZero();

    double numerator = 0;
    double denominator = 0;
    double r = 0;
    double exponent = 0;

    for ( int dim = 0; dim < valuesInterpolation.cols(); dim++ )
    {
      for ( int i = 0; i < m; i++ )
      {
        numerator = 0;
        denominator = 0;

        for ( int j = 0; j < n; j++ )
        {
          r = ( positionsInterpolation.row( i ) - positions.row( j ) ).norm();

          if ( r < 1.0e-13 )
            continue;

          exponent = std::pow( r, -c );
          numerator += values( j, dim ) * exponent;
          denominator += exponent;
        }

        valuesInterpolation( i, dim ) = numerator / denominator;
      }
    }
  }
}
