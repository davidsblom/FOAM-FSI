
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MultiLevelSolver_H
#define MultiLevelSolver_H

#include <iostream>
#include <memory>

#include "BaseMultiLevelSolver.H"
#include "RBFCoarsening.H"
#include "fvCFD.H"
#include "TPSFunction.H"

using std::shared_ptr;
using namespace rbf;

namespace fsi
{
    class MultiLevelSolver
    {
        public:
            MultiLevelSolver(
                shared_ptr<BaseMultiLevelSolver> solver,
                shared_ptr<BaseMultiLevelSolver> couplingGridSolver,
                int participantId,
                int level
                );

            MultiLevelSolver(
                shared_ptr<BaseMultiLevelSolver> solver,
                shared_ptr<BaseMultiLevelSolver> couplingGridSolver,
                shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh,
                shared_ptr<RBFCoarsening> rbfInterpToMesh,
                int participantId,
                int level
                );

            void interpToCouplingMesh(
                matrix & data,
                matrix & dataInterpolated
                );

            void interpToMesh(
                matrix & data,
                matrix & dataInterpolated
                );

            void solve(
                matrix & input,
                matrix & outputInterpolated
                );

            void updateInterpolationMatrices();

            shared_ptr<BaseMultiLevelSolver> solver;
            shared_ptr<BaseMultiLevelSolver> couplingGridSolver;
            shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
            shared_ptr<RBFCoarsening> rbfInterpToMesh;
            const int participantId;
            const int level;
            int couplingGridSize;
    };
}

#endif
