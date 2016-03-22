
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowSolidSolver_H
#define TubeFlowSolidSolver_H

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
    /*
     * The solid solver for the flexible tube fluid-solid interaction test case
     */
    class TubeFlowSolidSolver : public BaseMultiLevelSolver
    {
        public:
            TubeFlowSolidSolver(
                scalar a0,
                scalar cmk,
                scalar p0,
                scalar rho,
                scalar L,
                int N
                );

            virtual ~TubeFlowSolidSolver();

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

            void calcGrid();

            void solve(
                const fsi::vector & p,
                fsi::vector & a
                );

            scalar a0;
            scalar cmk;
            scalar p0;
            scalar rho;
            scalar dx;
            matrix grid;
    };
}


#endif
