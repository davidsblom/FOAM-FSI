
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef IDWInterpolation_H
#define IDWInterpolation_H

#include "BaseMultiLevelSolver.H"
#include <math.h>

namespace fsi
{
  class IDWInterpolation
  {
public:

    IDWInterpolation();

    IDWInterpolation ( double c );

    ~IDWInterpolation();

    void interpolate(
      const matrix & positions,
      const matrix & positionsInterpolation,
      const matrix & values,
      matrix & valuesInterpolation
      );

    double c;
  };
}

#endif
