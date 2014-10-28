
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MultiLevelFsiSolver_H
#define MultiLevelFsiSolver_H

#include "fvCFD.H"
#include "FsiSolver.H"
#include "MultiLevelSolver.H"

using namespace fsi;

class MultiLevelFsiSolver : public FsiSolver
{
public:

  MultiLevelFsiSolver(
    shared_ptr<MultiLevelSolver> fluidSolver,
    shared_ptr<MultiLevelSolver> solidSolver,
    shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures,
    bool parallel,
    int extrapolationOrder
    );

  ~MultiLevelFsiSolver();

  virtual void evaluate(
    const fsi::vector & input,
    fsi::vector & output,
    fsi::vector & R
    );

  void initTimeStep();

  bool isConvergence();

  bool isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    );

  void setSurrogateData(
    fsi::vector & xf,
    matrix & J
    );

  void setUseJacobian( bool useJacobian );

  shared_ptr<MultiLevelSolver> fluidSolver;
  shared_ptr<MultiLevelSolver> solidSolver;

  fsi::vector xf;
  matrix J;
  bool useJacobian;
  int iterCurrentTimeStep;
};

#endif
