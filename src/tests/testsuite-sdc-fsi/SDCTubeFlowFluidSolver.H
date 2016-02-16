
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

 #pragma once

 #include "TubeFlowFluidSolver.H"
 #include "SDCFsiSolverInterface.H"

 namespace tubeflow
 {
     class SDCTubeFlowFluidSolver : public TubeFlowFluidSolver, public sdc::SDCFsiSolverInterface
     {
     public:
         SDCTubeFlowFluidSolver(
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

        virtual ~SDCTubeFlowFluidSolver();

        virtual void evaluateFunction(
            const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            ) override;

        virtual void finalizeTimeStep() override;

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

        virtual void initTimeStep() override;

        virtual void setNumberOfImplicitStages( int k ) override;

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

        virtual void prepareImplicitSolve(
            bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs
        ) override;

        private:
            int k;

            std::deque<fsi::vector> pStages;
            std::deque<fsi::vector> uStages;
            std::deque<fsi::vector> aStages;
     };
 }
