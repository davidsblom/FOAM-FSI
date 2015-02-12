
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowLinearizedFluidSolver_H
#define TubeFlowLinearizedFluidSolver_H

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
  class TubeFlowLinearizedFluidSolver : public BaseMultiLevelSolver
  {
public:

    TubeFlowLinearizedFluidSolver(
      int N,
      double p0,
      double r0,
      double u0,
      double rho,
      double E,
      double h,
      double T,
      double dt,
      double L
      );

    virtual ~TubeFlowLinearizedFluidSolver();

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

    void factorizeMatrix();

    void solve(
      const fsi::vector & a,
      fsi::vector & p
      );

    void calcGrid();

    // Constants

    double p0;
    double r0;
    double u0;
    double a0;
    double rho;
    double h;
    double T;
    double dt;
    double dx;
    double tau;
    double cmk;
    double alpha;

    matrix grid;

    // Fields

    fsi::vector un;
    fsi::vector an;
    fsi::vector pn;
    fsi::vector u;
    fsi::vector a;
    fsi::vector p;

    Eigen::FullPivLU<matrix> lu;
  };
}

#endif
