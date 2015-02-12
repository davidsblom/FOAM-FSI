
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowLinearSolidSolver_H
#define TubeFlowLinearSolidSolver_H

#include "TubeFlowLinearizedSolidSolver.H"

using namespace fsi;

namespace tubeflow
{
  class TubeFlowLinearSolidSolver : public TubeFlowLinearizedSolidSolver
  {
public:

    TubeFlowLinearSolidSolver(
      int N,
      double nu,
      double rho,
      double h,
      double L,
      double dt,
      double G,
      double E0,
      double r0
      );

    virtual ~TubeFlowLinearSolidSolver();

    virtual void solve(
      const fsi::vector & p,
      fsi::vector & a
      );

    double h;
    double E0;
    double nu;
  };
}

#endif
