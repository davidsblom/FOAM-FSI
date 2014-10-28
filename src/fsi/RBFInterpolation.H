
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
      const std::string & function,
      const double radius
      );

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

    void interpolate(
      const matrix & positions,
      const matrix & positionsInterpolation,
      const matrix & values,
      matrix & valuesInterpolation
      );

    void interpolate(
      const matrix & positions,
      const matrix & positionsInterpolation,
      matrix & valuesInterpolation
      );

    double TPS( const double x );

    double WendlandC2( const double x );

    bool computed;
    int n_A;
    int n_B;
    int dimGrid;
    double radius;
    std::string function;
    matrix Hhat;
    matrix B;
  };
}

#endif
