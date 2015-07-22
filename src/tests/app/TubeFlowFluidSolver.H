
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TubeFlowFluidSolver_H
#define TubeFlowFluidSolver_H

#include "BaseMultiLevelSolver.H"
#include <unsupported/Eigen/NonLinearOptimization>

using namespace fsi;

namespace tubeflow
{
    class TubeFlowFluidSolver : public BaseMultiLevelSolver
    {
public:

        TubeFlowFluidSolver(
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

        virtual ~TubeFlowFluidSolver();

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

        bool isConvergence( const fsi::vector & R );

        void solve(
            const fsi::vector & input,
            fsi::vector & output
            );

        scalar evaluateInletVelocityBoundaryCondition();

        scalar evaluateOutputPressureBoundaryCondition(
            scalar pout_n,
            scalar uout_n,
            scalar uout
            );

        void evaluateJacobian(
            const fsi::vector & x,
            const fsi::vector & a,
            const fsi::vector & un,
            const fsi::vector & pn,
            const fsi::vector & an,
            matrix & J
            );

        void evaluateResidual(
            const fsi::vector & x,
            const fsi::vector & a,
            const fsi::vector & un,
            const fsi::vector & pn,
            const fsi::vector & an,
            fsi::vector & R
            );

        scalar a0;
        scalar u0;
        scalar p0;
        scalar dt;
        scalar dx;
        scalar cmk;
        scalar rho;
        scalar T;
        scalar alpha;
        scalar tau;
        scalar p_outn;
        scalar p_out;

        fsi::vector un;
        fsi::vector pn;
        fsi::vector an;
        fsi::vector u;
        fsi::vector p;
        fsi::vector a;

        int iter;
        int minIter;
        int maxIter;
        scalar tol;
        int nbRes;
        int nbJac;

        matrix grid;
    };
}

#endif
