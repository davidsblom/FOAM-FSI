
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include "BaseMultiLevelSolver.H"
#include <math.h>
#include <string>

namespace fsi
{
  class RBFInterpolation
  {
public:

    RBFInterpolation();

    RBFInterpolation(
      const matrix & positions,
      const matrix & positionsInterpolation,
      const std::string & function,
      const double radius
      );

    void compute(
      const matrix & positions,
      const matrix & positionsInterpolation
      );

    double evaluateFunction( const double x );

    void interpolate(
      const matrix & values,
      matrix & valuesInterpolation
      );

    double TPS( const double x );

    bool computed;
    int n_A;
    int n_B;
    int dimGrid;
    double radius;
    std::string function;
    matrix Hhat;
  };
}

#endif
