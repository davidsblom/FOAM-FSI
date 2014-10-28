
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MLIQNILSSolver_H
#define MLIQNILSSolver_H

#include "ImplicitMultiLevelFsiSolver.H"
#include "Solver.H"
#include <deque>

using namespace fsi;

class MLIQNILSSolver : public Solver
{
public:

  MLIQNILSSolver(
    shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models,
    bool synchronization
    );

  virtual void run();

  void finalizeTimeStep();

  void initTimeStep();

  bool isRunning();

  void solve();

  void solveTimeStep();

  shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;
  shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
  bool init;
  bool synchronization;
};

#endif
