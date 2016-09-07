
/*
 * Copyright [2016] <David Blom>
 */

#ifndef TubeFlowSolidSolver_H
#define TubeFlowSolidSolver_H

#include "BaseMultiLevelSolver.H"

namespace tubeflow
{
    /*
     * The solid solver for the flexible tube fluid-solid interaction test case
     */
    class TubeFlowSolidSolver : public fsi::BaseMultiLevelSolver
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

            virtual void getReadPositions( fsi::matrix & readPositions );

            virtual void getWritePositions( fsi::matrix & writePositions );

            virtual void initTimeStep();

            virtual bool isRunning();

            virtual void resetSolution();

            virtual void solve(
                const fsi::matrix & input,
                fsi::matrix & output
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
            fsi::matrix grid;
    };
} // namespace tubeflow


#endif
