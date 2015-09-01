
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowLinearSolidSolver_H
#define TubeFlowLinearSolidSolver_H

#include "TubeFlowLinearizedSolidSolver.H"

using namespace fsi;

namespace tubeflow
{
    class TubeFlowLinearSolidSolver : public TubeFlowLinearizedSolidSolver
    {
public:

        TubeFlowLinearSolidSolver(
            int N,
            scalar nu,
            scalar rho,
            scalar h,
            scalar L,
            scalar dt,
            scalar G,
            scalar E0,
            scalar r0
            );

        virtual ~TubeFlowLinearSolidSolver();

        virtual void solve(
            const fsi::vector & p,
            fsi::vector & a
            );

        scalar h;
        scalar E0;
        scalar nu;
    };
}

#endif
