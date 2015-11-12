
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
    class dealiiSolidSolver : public dealiifsi::LinearElasticity<dimension>, public BaseMultiLevelSolver
    {
public:

        dealiiSolidSolver( dealiifsi::DataStorage & data );

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
