
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ResidualRelativeConvergenceMeasure.H"

namespace fsi
{
    ResidualRelativeConvergenceMeasure::ResidualRelativeConvergenceMeasure(
        int dataId,
        bool suffices,
        scalar convergenceLimit
        )
        :
        ConvergenceMeasure( dataId, suffices ),
        isFirstIteration( true ),
        convergenceLimit( convergenceLimit ),
        normFirstResidual( std::numeric_limits<double>::max() ),
        normDiff( 0 )
    {
        assert( convergenceLimit > 0 );
        assert( convergenceLimit < 1 );
    }

    void ResidualRelativeConvergenceMeasure::measure(
        vector & oldValues,
        vector & newValues
        )
    {
        normDiff = (newValues - oldValues).norm();

        if ( isFirstIteration )
        {
            normFirstResidual = normDiff;
            isFirstIteration = false;
        }

        isConvergence_ = normDiff < normFirstResidual * convergenceLimit;
    }

    void ResidualRelativeConvergenceMeasure::newMeasurementSeries()
    {
        isConvergence_ = false;
        isFirstIteration = true;
        normFirstResidual = std::numeric_limits<double>::max();
    }

    bool ResidualRelativeConvergenceMeasure::isConvergence()
    {
        return isConvergence_;
    }

    void ResidualRelativeConvergenceMeasure::printState()
    {
        Info << "residual relative convergence measure: ";
        Info << "two-norm diff = " << normDiff;
        Info << ", limit = " << normFirstResidual * convergenceLimit;
        Info << ", suffices = ";

        if ( suffices() )
            Info << "true";
        else
            Info << "false";

        Info << ", conv = ";

        if ( isConvergence_ )
            Info << "true";
        else
            Info << "false";

        Info << endl;
    }
}
