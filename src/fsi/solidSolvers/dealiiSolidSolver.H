
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <LinearElasticity.h>
#include "BaseMultiLevelSolver.H"

namespace fsi
{
    template <int dimension>
    class dealiiSolidSolver : public Step23::LinearElasticity<dimension>, public BaseMultiLevelSolver
    {
public:

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

        virtual void finalizeTimeStep();

        virtual void getReadPositions( matrix & readPositions );

        virtual void getWritePositions( matrix & writePositions );

        virtual void initTimeStep();

        virtual bool isRunning();

        virtual void resetSolution();

        virtual void solve(
            const matrix & input,
            matrix & output
            );
    };

    #include "solidSolvers/dealiiSolidSolver.tpp"
}
