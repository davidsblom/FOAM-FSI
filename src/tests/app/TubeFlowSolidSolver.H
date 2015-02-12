
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowSolidSolver_H
#define TubeFlowSolidSolver_H

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
  /*
   * The solid solver for the flexible tube fluid-solid interaction test case
   */
  class TubeFlowSolidSolver : public BaseMultiLevelSolver
  {
public:

    TubeFlowSolidSolver(
      double a0,
      double cmk,
      double p0,
      double rho,
      double L,
      int N
      );

    virtual ~TubeFlowSolidSolver();

    virtual void finalizeTimeStep();

    virtual void getReadPositions( matrix & readPositions );

    virtual void getWritePositions( matrix & writePositions );

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void resetSolution();

    virtual void solve(
      const matrix & input,
      matrix & output
      );

    void calcGrid();

    void solve(
      const fsi::vector & p,
      fsi::vector & a
      );

    double a0;
    double cmk;
    double p0;
    double rho;
    double dx;
    matrix grid;
  };
}


#endif
