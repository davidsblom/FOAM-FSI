
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "TubeFlowSolidSolver.H"
#include "SDCFsiSolverInterface.H"

namespace tubeflow {
using namespace fsi;

/*
 * The solid solver for the flexible tube fluid-solid interaction test case
 */
class SDCTubeFlowSolidSolver : public TubeFlowSolidSolver, public sdc::SDCFsiSolverInterface {
    public:
        SDCTubeFlowSolidSolver(scalar a0,
            scalar cmk,
            scalar p0,
            scalar rho,
            scalar L,
            int N
            );

        virtual ~SDCTubeFlowSolidSolver();

        virtual void evaluateFunction(const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            ) override;

        virtual void finalizeTimeStep() override;

        virtual int getDOF() override;

        virtual void getSolution(fsi::vector & solution,
            fsi::vector & f
            ) override;

        virtual void setSolution(const fsi::vector & solution,
            const fsi::vector & f
            ) override;

        virtual scalar getEndTime() override;

        virtual scalar getTimeStep() override;

        virtual void nextTimeStep() override;

        virtual void initTimeStep() override;

        virtual void setNumberOfImplicitStages(int k) override;

        virtual void implicitSolve(bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs,
            fsi::vector & f,
            fsi::vector & result
            ) override;

        virtual void prepareImplicitSolve(bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs
            ) override;

        virtual void getVariablesInfo(std::deque<int> & dof,
            std::deque<bool> & enabled,
            std::deque<std::string> & names
            ) override;
};
}
