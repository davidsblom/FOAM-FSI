
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "TubeFlowLinearizedSolidSolver.H"

namespace tubeflow
{
    using namespace fsi;

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
                scalar r0,
                scalar T
                );

            virtual ~TubeFlowLinearSolidSolver();

            virtual void solve(
                const fsi::vector & p,
                fsi::vector & a
                );

            const scalar E0;
            const scalar nu;
            const scalar T;
    };
}
