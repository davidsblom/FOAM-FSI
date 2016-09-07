
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <unsupported/Eigen/NonLinearOptimization>
#include "BaseMultiLevelSolver.H"

namespace tubeflow {
class TubeFlowFluidSolver : public fsi::BaseMultiLevelSolver {
    public:
        TubeFlowFluidSolver(scalar a0,
            scalar u0,
            scalar p0,
            scalar dt,
            scalar cmk,
            int N,
            scalar L,
            scalar T,
            scalar rho
            );

        TubeFlowFluidSolver(scalar a0,
            scalar u0,
            scalar p0,
            scalar dt,
            scalar cmk,
            int N,
            scalar L,
            scalar T,
            scalar rho,
            bool diffJacobian
            );

        TubeFlowFluidSolver(scalar a0,
            scalar u0,
            scalar p0,
            scalar dt,
            scalar cmk,
            int N,
            scalar L,
            scalar T,
            scalar rho,
            scalar tol
            );

        virtual ~TubeFlowFluidSolver();

        virtual void finalizeTimeStep();

        virtual void getReadPositions(fsi::matrix & readPositions);

        virtual void getWritePositions(fsi::matrix & writePositions);

        virtual void initTimeStep();

        virtual bool isRunning();

        virtual void resetSolution();

        virtual void solve(const fsi::matrix & input,
            fsi::matrix & output
            );

        void calcGrid();

        bool isConvergence(const fsi::vector & R);

        void solve(const fsi::vector & input,
            fsi::vector & output
            );

        virtual scalar evaluateInletVelocityBoundaryCondition();

        virtual scalar evaluateOutputPressureBoundaryCondition(const fsi::vector & u);

        void evaluateJacobian(const fsi::vector & x,
            const fsi::vector & a,
            const fsi::vector & un,
            const fsi::vector & pn,
            const fsi::vector & an,
            fsi::matrix & J
            );

        void evaluateResidual(const fsi::vector & x,
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
        scalar L;
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
        fsi::vector rhs;

        int iter;
        int minIter;
        int maxIter;
        scalar tol;
        int nbRes;
        int nbJac;

        fsi::matrix grid;

        bool diffJacobian;
};
} // namespace tubeflow
