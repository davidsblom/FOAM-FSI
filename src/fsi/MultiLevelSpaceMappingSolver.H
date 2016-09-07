
/*
 * Copyright [2016] <David Blom>
 */

#ifndef MultiLevelSpaceMappingSolver_H
#define MultiLevelSpaceMappingSolver_H

#include <deque>
#include <memory>
#include "Solver.H"
#include "SpaceMappingSolver.H"

namespace fsi {
class MultiLevelSpaceMappingSolver : public Solver {
    public:
        MultiLevelSpaceMappingSolver(shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > solvers,
            shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models,
            bool synchronization
            );

        virtual void run();

        void initTimeStep();

        void finalizeTimeStep();

        bool isRunning();

        void solve();

        void solveTimeStep();

        shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > solvers;
        shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;
        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
        bool init;
        const bool synchronization;
};
} // namespace fsi

#endif
