
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFCoarsening_H
#define RBFCoarsening_H

#include "RBFInterpolation.H"
#include "fvCFD.H"

namespace rbf
{
  class RBFCoarsening
  {
public:

    RBFCoarsening(
      std::shared_ptr<RBFInterpolation> rbf,
      bool enabled,
      double tol,
      int coarseningMinPoints,
      int coarseningMaxPoints
      );

    void compute(
      const matrix & positions,
      const matrix & positionsInterpolation
      );

    void interpolate(
      const matrix & values,
      matrix & valuesInterpolation
      );

    void setNbMovingAndStaticFaceCenters(
      int nbMovingFaceCenters,
      int nbStaticFaceCenters
      );

    std::shared_ptr<RBFInterpolation> rbf;
    bool enabled;
    double tol;
    int coarseningMinPoints;
    int coarseningMaxPoints;
    Eigen::VectorXi selectedPositions;
    int nbStaticFaceCentersRemove;
  };
}

#endif
