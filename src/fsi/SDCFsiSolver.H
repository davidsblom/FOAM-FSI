
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCFsiSolver_H
#define SDCFsiSolver_H

#include "SDCFsiSolverInterface.H"
#include "PostProcessing.H"

using namespace sdc;

namespace fsi
{
    class SDCFsiSolver : public SDCSolver
    {
public:

        SDCFsiSolver(
            std::shared_ptr<SDCFsiSolverInterface> fluid,
            std::shared_ptr<SDCFsiSolverInterface> solid,
            std::shared_ptr<PostProcessing> postProcessing
            );

        virtual ~SDCFsiSolver();

        virtual void evaluateFunction(
            const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            );

        virtual void finalizeTimeStep();

        virtual int getDOF();

        virtual scalar getScalingFactor();

        virtual void getSolution(
            fsi::vector & solution,
            fsi::vector & f
            );

        virtual void setSolution(
            const fsi::vector & solution,
            const fsi::vector & f
            );

        virtual scalar getEndTime();

        virtual scalar getTimeStep();

        virtual void nextTimeStep();

        virtual void initTimeStep();

        virtual void setNumberOfImplicitStages( int k );

        virtual void implicitSolve(
            bool corrector,
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

        virtual void getVariablesInfo(
            std::deque<int> & dof,
            std::deque<bool> & enabled,
            std::deque<std::string> & names
            );

        virtual bool isConverged();

private:

        std::shared_ptr<SDCFsiSolverInterface> fluid;
        std::shared_ptr<SDCFsiSolverInterface> solid;
        std::shared_ptr<PostProcessing> postProcessing;
        int dofFluid;
        int dofSolid;
    };
}

#endif
