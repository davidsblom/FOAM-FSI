
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
    class TubeFlowExplicitLinearSolidSolver : public BaseMultiLevelSolver
    {
        public:
            TubeFlowExplicitLinearSolidSolver(
                int N,
                scalar nu,
                scalar rho,
                scalar h,
                scalar L,
                scalar dt,
                scalar G,
                scalar E0,
                scalar r0,
                scalar p0,
                scalar T
                );

            virtual ~TubeFlowExplicitLinearSolidSolver();

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

            virtual void solve(
                const fsi::vector & p,
                fsi::vector & a
                );

            void calcGrid();

            fsi::vector evaluateFunction(
                const fsi::vector & rn,
                const fsi::vector & pn,
                const fsi::vector & un
                );

            // Constants

            const int N;
            scalar dt;
            const scalar r0;
            const scalar kappa;
            const scalar dx;
            const scalar G;
            const scalar E0;
            const scalar nu;
            const scalar h;
            const scalar rho;
            const scalar T;

            matrix grid;

            // Fields

            fsi::vector un;
            fsi::vector rn;
            fsi::vector u;
            fsi::vector r;
            fsi::vector rhs;
            fsi::vector p;
            fsi::vector pn;
    };
}
