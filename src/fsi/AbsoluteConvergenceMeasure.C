
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
        isConvergence = normDiff < convergenceLimit;
    }

    void AbsoluteConvergenceMeasure::newMeasurementSeries()
    {
        isConvergence = false;
    }

    bool AbsoluteConvergenceMeasure::isConvergence()
    {
        return isConvergence;
    }

    void AbsoluteConvergenceMeasure::printState()
    {
        Info << "absolute convergence measure: ";
        Info << "two-norm diff = " << normDiff;
        Info << ", limit = " << convergenceLimit;
        Info << ", conv = ";
        if ( isConvergence )
            Info << "true";
        else
            Info << "false";
        Info << endl;
    }

}
