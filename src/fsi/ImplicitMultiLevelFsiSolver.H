
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ImplicitMultiLevelFsiSolver_H
#define ImplicitMultiLevelFsiSolver_H

#include "PostProcessing.H"
#include "MultiLevelFsiSolver.H"
#include "Solver.H"
#include "SurrogateModel.H"

using namespace fsi;

class ImplicitMultiLevelFsiSolver : public SurrogateModel, public Solver
{
public:

  ImplicitMultiLevelFsiSolver(
    shared_ptr<MultiLevelFsiSolver> fsi,
    shared_ptr<PostProcessing> postProcessing
    );

  virtual ~ImplicitMultiLevelFsiSolver();

  virtual bool allConverged();

  virtual void evaluate(
    const fsi::vector & x,
    fsi::vector & output,
    fsi::vector & R
    );

  virtual bool isConvergence();

  virtual bool isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    );

  virtual void optimize(
    const fsi::vector & x0,
    fsi::vector & xk
    );

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk
    );

  virtual void run();

  virtual void setSurrogateData(
    fsi::vector & xf,
    matrix & J
    );

  virtual void setUseJacobian( bool useJacobian );

  void finalizeTimeStep();

  void initTimeStep();

  void solve();

  void solveTimeStep();

  shared_ptr<MultiLevelFsiSolver> fsi;
  shared_ptr<PostProcessing> postProcessing;
  bool init;
};

#endif
