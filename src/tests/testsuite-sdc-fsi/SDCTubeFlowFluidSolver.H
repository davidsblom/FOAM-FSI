
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

 #pragma once

 #include "TubeFlowFluidSolver.H"
 #include "SDCSolver.H"

 namespace tubeflow
 {
     class SDCTubeFlowFluidSolver : public TubeFlowFluidSolver, public sdc::SDCSolver
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
            );

        virtual void finalizeTimeStep();

        virtual int getDOF();

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

        private:
            int k;

            std::deque<fsi::vector> pStages;
            std::deque<fsi::vector> uStages;
            std::deque<fsi::vector> aStages;
     };
 }
