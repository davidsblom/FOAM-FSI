
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RelativeConvergenceMeasure.H"

namespace fsi
{
    RelativeConvergenceMeasure::RelativeConvergenceMeasure(
        int dataId,
        double convergenceLimit
        )
        :
        ConvergenceMeasure( dataId, convergenceLimit ),
        normDiff( 0 ),
        norm( 0 ),
        convergenceLimit( convergenceLimit )
    {
        assert( convergenceLimit > 0 );
        assert( convergenceLimit < 1 );

        epsilon = std::sqrt( SMALL );
    }

    void RelativeConvergenceMeasure::measure(
        vector & oldValues,
        vector & newValues
        )
    {
        assert( oldValues.rows() == newValues.rows() );
        assert( oldValues.rows() > 0 );

        vector diff = newValues - oldValues;
        normDiff = diff.norm();
        norm = newValues.norm();

        isConvergence_ = (epsilon + norm) * convergenceLimit > normDiff;
    }

    void RelativeConvergenceMeasure::newMeasurementSeries()
    {
        isConvergence_ = false;
    }

    bool RelativeConvergenceMeasure::isConvergence()
    {
        return isConvergence_;
    }

    void RelativeConvergenceMeasure::printState()
    {
        Info << "relative convergence measure: two-norm diff = "
             << normDiff
             << ", limit = "
             << (epsilon + norm) * convergenceLimit
             << ", conv = ";

        if ( isConvergence_ )
            Info << "true";

        if ( !isConvergence_ )
            Info << "false";

        Info << endl;
    }
}
