
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MonolithicFsiSolver_H
#define MonolithicFsiSolver_H

#include "BaseMultiLevelSolver.H"
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace fsi;

namespace tubeflow
{
  class MonolithicFsiSolver
  {
public:

    MonolithicFsiSolver(
      double a0,
      double u0,
      double p0,
      double dt,
      double cmk,
      int N,
      double L,
      double T,
      double rho
      );

    double evaluateInletVelocityBoundaryCondition();

    double evaluateOutputPressureBoundaryCondition(
      double pout_n,
      double uout_n,
      double uout
      );

    void evaluateResidual(
      fsi::vector & x,
      fsi::vector & un,
      fsi::vector & pn,
      fsi::vector & an,
      fsi::vector & R
      );

    void finalizeTimeStep();

    void initTimeStep();

    bool isRunning();

    void run();

    void solve();

    void solveTimeStep();

    double a0;
    double u0;
    double p0;
    double dt;
    double cmk;
    int N;
    double L;
    double T;
    double rho;

    fsi::vector un;
    fsi::vector pn;
    fsi::vector an;
    fsi::vector u;
    fsi::vector p;
    fsi::vector a;

    bool init;

    double dx;
    double alpha;
    double tau;

    int timeIndex;
    double t;

    double p_outn;
    double p_out;
  };
}

#endif
