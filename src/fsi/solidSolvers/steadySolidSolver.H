
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef steadySolidSolver_H
#define steadySolidSolver_H

#include "SolidSolver.H"

using namespace fsi;

class steadySolidSolver : public SolidSolver
{
public:

  steadySolidSolver(
    string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    );

  virtual ~steadySolidSolver();

  virtual void solve();
};

#endif
