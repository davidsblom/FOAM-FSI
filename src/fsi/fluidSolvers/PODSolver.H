
/*
 * Class
 *   FluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef PODSolver_H
#define PODSolver_H

#include "BaseMultiLevelSolver.H"
#include "RBFInterpolation.H"
#include "ImplicitMultiLevelFsiSolver.H"

using namespace fsi;
using std::shared_ptr;
using std::deque;

class PODSolver : public BaseMultiLevelSolver
{
public:

  PODSolver(
    shared_ptr<ImplicitMultiLevelFsiSolver> fsi,
    int participantId
    );

  virtual ~PODSolver();

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

  deque<fsi::vector> inputList;
  deque<fsi::vector> outputList;

private:

  shared_ptr<ImplicitMultiLevelFsiSolver> fsi;

  int participantId;
};

#endif
