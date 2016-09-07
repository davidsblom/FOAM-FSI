
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <memory>
#include <string>
#include <deque>
#include "SDCFsiSolverInterface.H"
#include "PostProcessing.H"

namespace fsi {
class SDCFsiSolver : public sdc::SDCSolver {
    public:
        SDCFsiSolver(std::shared_ptr<sdc::SDCFsiSolverInterface> fluid,
            std::shared_ptr<sdc::SDCFsiSolverInterface> solid,
            std::shared_ptr<PostProcessing> postProcessing,
            int extrapolationOrder
            );

        virtual ~SDCFsiSolver();

        virtual void evaluateFunction(const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            );

        virtual void finalizeTimeStep();

        virtual int getDOF();

        virtual void getSolution(fsi::vector & solution,
            fsi::vector & f
            );

        virtual void setSolution(const fsi::vector & solution,
            const fsi::vector & f
            );

        virtual scalar getEndTime();

        virtual scalar getTimeStep();

        virtual void nextTimeStep();

        virtual void initTimeStep();

        virtual void setNumberOfImplicitStages(int k);

        virtual void implicitSolve(bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs,
            fsi::vector & f,
            fsi::vector & result
            );

        virtual scalar getStartTime();

        virtual void getVariablesInfo(std::deque<int> & dof,
            std::deque<bool> & enabled,
            std::deque<std::string> & names
            );

        virtual bool isConverged();

    private:
        std::shared_ptr<sdc::SDCFsiSolverInterface> fluid;
        std::shared_ptr<sdc::SDCFsiSolverInterface> solid;
        std::shared_ptr<PostProcessing> postProcessing;
        int dofFluid;
        int dofSolid;

        int k;
        std::deque<fsi::vector> xStages;
        std::deque<fsi::vector> xStagesPrevious;
        int extrapolationOrder;
};
} // namespace fsi
