
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowLinearizedSolidSolver_H
#define TubeFlowLinearizedSolidSolver_H

#include "BaseMultiLevelSolver.H"

using namespace fsi;

namespace tubeflow
{
    class TubeFlowLinearizedSolidSolver : public BaseMultiLevelSolver
    {
public:

        TubeFlowLinearizedSolidSolver(
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

        virtual ~TubeFlowLinearizedSolidSolver();

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

        void factorizeMatrix();

        void calcGrid();

        // Constants

        int N;
        scalar dt;
        scalar r0;
        scalar kappa;
        scalar dx;
        scalar alpha;
        scalar beta;
        scalar gamma;

        matrix grid;

        // Fields

        fsi::vector un;
        fsi::vector rn;
        fsi::vector u;
        fsi::vector r;

        Eigen::FullPivLU<matrix> lu;
    };
}

#endif
