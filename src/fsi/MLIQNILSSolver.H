
/*
 * Copyright [2016] <David Blom>
 */

#ifndef MLIQNILSSolver_H
#define MLIQNILSSolver_H

#include "ImplicitMultiLevelFsiSolver.H"
#include "Solver.H"
#include <deque>

namespace fsi {

class MLIQNILSSolver : public Solver
{
    public:
        MLIQNILSSolver(
            shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models,
            bool synchronization
            );

        virtual void run();

        void finalizeTimeStep();

        void initTimeStep();

        bool isRunning();

        void solve();

        void solveTimeStep();

        shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;
        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
        bool init;
        bool synchronization;
};

} // namespace fsi

#endif
