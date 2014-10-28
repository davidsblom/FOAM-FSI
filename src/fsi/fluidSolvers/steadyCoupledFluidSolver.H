
/*
 * Class
 *   steadyCoupledFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef steadyCoupledFluidSolver_H
#define steadyCoupledFluidSolver_H

#include "CoupledFluidSolver.H"

using namespace fsi;

class steadyCoupledFluidSolver : public CoupledFluidSolver
{
public:

  steadyCoupledFluidSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    );

  virtual ~steadyCoupledFluidSolver();

  virtual void solve();
};

#endif
