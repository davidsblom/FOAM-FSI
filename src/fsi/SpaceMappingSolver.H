
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SpaceMappingSolver_H
#define SpaceMappingSolver_H

#include "ImplicitMultiLevelFsiSolver.H"
#include "SpaceMapping.H"

using namespace fsi;

class SpaceMappingSolver : public ImplicitMultiLevelFsiSolver
{
public:

  SpaceMappingSolver(
    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel,
    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel,
    shared_ptr<SpaceMapping> spaceMapping
    );

  virtual ~SpaceMappingSolver();

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

  bool isRunning();

  void solve();

  void solveTimeStep();

  shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
  shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel;
  shared_ptr<SpaceMapping> spaceMapping;
  bool init;
};

#endif
