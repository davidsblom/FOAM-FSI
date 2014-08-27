
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MultiLevelSolver_H
#define MultiLevelSolver_H

#include <iostream>
#include <memory>

#include "BaseMultiLevelSolver.H"
#include "RBFInterpolation.H"
#include "fvCFD.H"

using std::shared_ptr;

namespace fsi
{
  class MultiLevelSolver
  {
public:

    MultiLevelSolver(
      shared_ptr<BaseMultiLevelSolver> solver,
      shared_ptr<BaseMultiLevelSolver> couplingGridSolver,
      int participantId,
      int level
      );

    void interpToCouplingMesh(
      matrix & data,
      matrix & dataInterpolated
      );

    void interpToMesh(
      matrix & data,
      matrix & dataInterpolated
      );

    void solve(
      matrix & input,
      matrix & outputInterpolated
      );

    void updateInterpolationMatrices();

    shared_ptr<BaseMultiLevelSolver> solver;
    shared_ptr<BaseMultiLevelSolver> couplingGridSolver;
    shared_ptr<RBFInterpolation> rbfInterpToCouplingMesh;
    shared_ptr<RBFInterpolation> rbfInterpToMesh;
    int participantId;
    int level;
    int couplingGridSize;
  };
}

#endif
