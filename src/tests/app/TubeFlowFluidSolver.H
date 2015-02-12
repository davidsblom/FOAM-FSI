
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowFluidSolver_H
#define TubeFlowFluidSolver_H

#include "BaseMultiLevelSolver.H"
#include <unsupported/Eigen/NonLinearOptimization>

using namespace fsi;

namespace tubeflow
{
  class TubeFlowFluidSolver : public BaseMultiLevelSolver
  {
public:

    TubeFlowFluidSolver(
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

    virtual ~TubeFlowFluidSolver();

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

    bool isConvergence( const fsi::vector & R );

    void solve(
      const fsi::vector & input,
      fsi::vector & output
      );

    double evaluateInletVelocityBoundaryCondition();

    double evaluateOutputPressureBoundaryCondition(
      double pout_n,
      double uout_n,
      double uout
      );

    void evaluateJacobian(
      const fsi::vector & x,
      const fsi::vector & a,
      const fsi::vector & un,
      const fsi::vector & pn,
      const fsi::vector & an,
      matrix & J
      );

    void evaluateResidual(
      const fsi::vector & x,
      const fsi::vector & a,
      const fsi::vector & un,
      const fsi::vector & pn,
      const fsi::vector & an,
      fsi::vector & R
      );

    double a0;
    double u0;
    double p0;
    double dt;
    double dx;
    double cmk;
    double rho;
    double T;
    double alpha;
    double tau;
    double p_outn;
    double p_out;

    fsi::vector un;
    fsi::vector pn;
    fsi::vector an;
    fsi::vector u;
    fsi::vector p;
    fsi::vector a;

    int iter;
    int minIter;
    int maxIter;
    double tol;
    int nbRes;
    int nbJac;

    matrix grid;
  };
}

#endif
