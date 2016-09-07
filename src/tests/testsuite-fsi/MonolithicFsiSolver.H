
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

namespace tubeflow {
class MonolithicFsiSolver {
    public:
        MonolithicFsiSolver(const scalar a0,
            const scalar u0,
            const scalar p0,
            const scalar dt,
            const scalar cmk,
            const int N,
            const scalar L,
            const scalar T,
            const scalar rho
            );

        scalar evaluateInletVelocityBoundaryCondition();

        scalar evaluateOutputPressureBoundaryCondition(const scalar pout_n,
            const scalar uout_n,
            const scalar uout
            );

        void evaluateResidual(const fsi::vector & x,
            const fsi::vector & un,
            const fsi::vector & pn,
            const fsi::vector & an,
            fsi::vector & R
            );

        void finalizeTimeStep();

        void initTimeStep();

        bool isRunning() const;

        void run();

        void solve();

        void solveTimeStep();

        const scalar a0;
        const scalar u0;
        const scalar p0;
        const scalar dt;
        const scalar cmk;
        const int N;
        const scalar L;
        const scalar T;
        const scalar rho;

        fsi::vector un;
        fsi::vector pn;
        fsi::vector an;
        fsi::vector u;
        fsi::vector p;
        fsi::vector a;

        bool init;

        const scalar dx;
        const scalar alpha;
        const scalar tau;

        int timeIndex;
        scalar t;

        scalar p_outn;
        scalar p_out;
};
}

#endif
