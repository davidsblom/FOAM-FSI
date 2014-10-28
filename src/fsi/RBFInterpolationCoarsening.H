
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFInterpolationCoarsening_H
#define RBFInterpolationCoarsening_H

#include "RBFInterpolation.H"
#include <deque>

namespace fsi
{
  class RBFInterpolationCoarsening
  {
public:

    RBFInterpolationCoarsening( double tol );

    ~RBFInterpolationCoarsening();

    void coarsen(
      const matrix & positions,
      const matrix & values,
      matrix & selectedPositions
      );

    void interpolate(
      const matrix & positions,
      const matrix & positionsInterpolation,
      const matrix & values,
      matrix & valuesInterpolation
      );

    double tol;

    RBFInterpolation rbf;
    matrix positions;
    matrix positionsInterpolation;
  };
}

#endif
