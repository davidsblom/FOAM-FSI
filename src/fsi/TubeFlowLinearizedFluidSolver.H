
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowLinearizedFluidSolver_H
#define TubeFlowLinearizedFluidSolver_H

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
    class TubeFlowLinearizedFluidSolver : public BaseMultiLevelSolver
    {
        public:
            TubeFlowLinearizedFluidSolver(
                int N,
                scalar p0,
                scalar r0,
                scalar u0,
                scalar rho,
                scalar E,
                scalar h,
                scalar T,
                scalar dt,
                scalar L
                );

            virtual ~TubeFlowLinearizedFluidSolver();

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

            void factorizeMatrix();

            void solve(
                const fsi::vector & a,
                fsi::vector & p
                );

            void calcGrid();

            // Constants

            scalar p0;
            scalar r0;
            scalar u0;
            scalar a0;
            scalar rho;
            scalar h;
            scalar T;
            scalar dt;
            scalar dx;
            scalar tau;
            scalar cmk;
            scalar alpha;

            matrix grid;

            // Fields

            fsi::vector un;
            fsi::vector an;
            fsi::vector pn;
            fsi::vector u;
            fsi::vector a;
            fsi::vector p;

            Eigen::FullPivLU<matrix> lu;
    };
}

#endif
