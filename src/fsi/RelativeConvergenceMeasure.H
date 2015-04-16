
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RelativeConvergenceMeasure_H
#define RelativeConvergenceMeasure_H

#include "ConvergenceMeasure.H"
#include "fvCFD.H"

namespace fsi
{
    class RelativeConvergenceMeasure : public ConvergenceMeasure
    {
public:

        RelativeConvergenceMeasure(
            int dataId,
            double convergenceLimit
            );

        virtual void measure(
            vector & oldValues,
            vector & newValues
            );

        virtual void newMeasurementSeries();

        virtual bool isConvergence();

        virtual void printState();

        double normDiff;
        double norm;
        double convergenceLimit;
        double epsilon;
    };
}

#endif
