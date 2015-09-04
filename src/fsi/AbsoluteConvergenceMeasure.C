
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AbsoluteConvergenceMeasure.H"

namespace fsi
{
    AbsoluteConvergenceMeasure::AbsoluteConvergenceMeasure(
        int dataId,
        scalar convergenceLimit
        )
        :
        ConvergenceMeasure( dataId ),
        convergenceLimit( convergenceLimit ),
        normDiff( 0 )
    {
        assert( convergenceLimit > 0 );
        assert( convergenceLimit < 1 );
    }

    void AbsoluteConvergenceMeasure::measure(
        vector & oldValues,
        vector & newValues
        )
    {
        normDiff = (newValues - oldValues).norm();
        isConvergence_ = normDiff < convergenceLimit;
    }

    void AbsoluteConvergenceMeasure::newMeasurementSeries()
    {
        isConvergence_ = false;
    }

    bool AbsoluteConvergenceMeasure::isConvergence()
    {
        return isConvergence_;
    }

    void AbsoluteConvergenceMeasure::printState()
    {
        Info << "absolute convergence measure: ";
        Info << "two-norm diff = " << normDiff;
        Info << ", limit = " << convergenceLimit;
        Info << ", conv = ";
        if ( isConvergence_ )
            Info << "true";
        else
            Info << "false";
        Info << endl;
    }

}
