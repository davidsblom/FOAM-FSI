
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ResidualRelativeConvergenceMeasure_H
#define ResidualRelativeConvergenceMeasure_H

#include "ConvergenceMeasure.H"
#include "fvCFD.H"

namespace fsi
{
    class ResidualRelativeConvergenceMeasure : public ConvergenceMeasure
    {
public:

        ResidualRelativeConvergenceMeasure(
            int dataId,
            scalar convergenceLimit
            );

        virtual void measure(
            vector & oldValues,
            vector & newValues
            );

        virtual void newMeasurementSeries();

        virtual bool isConvergence();

        virtual void printState();

        bool isFirstIteration;
        scalar convergenceLimit;
        scalar normFirstResidual;
        scalar normDiff;
    };
}

#endif
