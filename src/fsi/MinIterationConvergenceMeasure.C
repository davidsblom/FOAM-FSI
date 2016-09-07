
/*
 * Copyright [2016] <David Blom>
 */

#include "fvCFD.H"
#include "MinIterationConvergenceMeasure.H"

namespace fsi
{
    MinIterationConvergenceMeasure::MinIterationConvergenceMeasure(
        int dataId,
        bool suffices,
        int minimumIterationCount
        )
        :
        ConvergenceMeasure( dataId, suffices ),
        minimumIterationCount( minimumIterationCount ),
        currentIteration( 0 )
    {
        assert( minimumIterationCount >= 0 );
    }

    void MinIterationConvergenceMeasure::measure(
        vector & /*oldValues*/,
        vector & /*newValues*/
        )
    {
        currentIteration++;

        isConvergence_ = minimumIterationCount <= currentIteration;
    }

    void MinIterationConvergenceMeasure::newMeasurementSeries()
    {
        currentIteration = 0;
        isConvergence_ = false;
    }

    bool MinIterationConvergenceMeasure::isConvergence()
    {
        return isConvergence_;
    }

    void MinIterationConvergenceMeasure::printState()
    {
        Info << "min iteration convergence measure: #it = "
             << currentIteration
             << " of "
             << minimumIterationCount
             << ", conv = ";

        if ( isConvergence_ )
            Info << "true";

        if ( !isConvergence_ )
            Info << "false";

        Info << endl;
    }
} // namespace fsi
