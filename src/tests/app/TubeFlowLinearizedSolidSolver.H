
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowLinearizedSolidSolver_H
#define TubeFlowLinearizedSolidSolver_H

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
  class TubeFlowLinearizedSolidSolver : public BaseMultiLevelSolver
  {
public:

    TubeFlowLinearizedSolidSolver(
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

    virtual ~TubeFlowLinearizedSolidSolver();

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

    virtual void solve(
      const fsi::vector & p,
      fsi::vector & a
      );

    void factorizeMatrix();

    void calcGrid();

    // Constants

    int N;
    double dt;
    double r0;
    double kappa;
    double dx;
    double alpha;
    double beta;
    double gamma;

    matrix grid;

    // Fields

    fsi::vector un;
    fsi::vector rn;
    fsi::vector u;
    fsi::vector r;

    Eigen::FullPivLU<matrix> lu;
  };
}

#endif
