
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MonolithicFsiSolver_H
#define MonolithicFsiSolver_H

#include "BaseMultiLevelSolver.H"
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace fsi;

namespace tubeflow
{
    class MonolithicFsiSolver
    {
public:

        MonolithicFsiSolver(
            scalar a0,
            scalar u0,
            scalar p0,
            scalar dt,
            scalar cmk,
            int N,
            scalar L,
            scalar T,
            scalar rho
            );

        scalar evaluateInletVelocityBoundaryCondition();

        scalar evaluateOutputPressureBoundaryCondition(
            scalar pout_n,
            scalar uout_n,
            scalar uout
            );

        void evaluateResidual(
            fsi::vector & x,
            fsi::vector & un,
            fsi::vector & pn,
            fsi::vector & an,
            fsi::vector & R
            );

        void finalizeTimeStep();

        void initTimeStep();

        bool isRunning();

        void run();

        void solve();

        void solveTimeStep();

        scalar a0;
        scalar u0;
        scalar p0;
        scalar dt;
        scalar cmk;
        int N;
        scalar L;
        scalar T;
        scalar rho;

        fsi::vector un;
        fsi::vector pn;
        fsi::vector an;
        fsi::vector u;
        fsi::vector p;
        fsi::vector a;

        bool init;

        scalar dx;
        scalar alpha;
        scalar tau;

        int timeIndex;
        scalar t;

        scalar p_outn;
        scalar p_out;
    };
}

#endif
