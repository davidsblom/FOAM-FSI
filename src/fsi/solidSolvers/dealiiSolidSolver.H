
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <mpi.h>
#include <LinearElasticity.h>
#include "BaseMultiLevelSolver.H"
#include "SDCFsiSolverInterface.H"

namespace fsi
{
    void copy(
        const dealii::Vector<double> & source,
        fsi::vector & target,
        unsigned int targetOffset
        );

    void copy(
        const fsi::vector & source,
        dealii::Vector<double> & target,
        unsigned int sourceOffset
        );

    template <int dimension>
    class dealiiSolidSolver : public dealiifsi::LinearElasticity<dimension>, public BaseMultiLevelSolver, public sdc::SDCFsiSolverInterface
    {
public:

        explicit dealiiSolidSolver( dealiifsi::DataStorage & data );

        dealiiSolidSolver(
            double time_step,
            double final_time,
            double theta,
            double degree,
            double gravity,
            double distributed_load,
            double rho,
            double E,
            double nu,
            unsigned int n_global_refines
            );

        virtual ~dealiiSolidSolver();

        virtual void finalizeTimeStep() override;

        virtual void getReadPositions( matrix & readPositions ) override;

        virtual void getWritePositions( matrix & writePositions ) override;

        virtual void initTimeStep() override;

        virtual bool isRunning() override;

        virtual void resetSolution() override;

        virtual void solve(
            const matrix & input,
            matrix & output
            ) override;

        virtual void evaluateFunction(
            const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            ) override;

        virtual int getDOF() override;

        virtual void getSolution(
            fsi::vector & solution,
            fsi::vector & f
            ) override;

        virtual void setSolution(
            const fsi::vector & solution,
            const fsi::vector & f
            ) override;

        virtual scalar getEndTime() override;

        virtual scalar getTimeStep() override;

        virtual void nextTimeStep() override;

        virtual void setNumberOfImplicitStages( int k ) override;

        virtual void prepareImplicitSolve(
            bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs
            ) override;

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
            ) override;

        virtual void finalizeImplicitSolve( int k ) override;

        virtual void getVariablesInfo(
            std::deque<int> & dof,
            std::deque<bool> & enabled,
            std::deque<std::string> & names
            ) override;

private:

        int k, kindex;
        std::deque<dealii::Vector<double> > UStages;
        std::deque<dealii::Vector<double> > VStages;
    };

    #include "solidSolvers/dealiiSolidSolver.tpp"
}
