
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef MinIterationConvergenceMeasure_H
#define MinIterationConvergenceMeasure_H

#include "ConvergenceMeasure.H"

namespace fsi
{
    class MinIterationConvergenceMeasure : public ConvergenceMeasure
    {
        public:
            MinIterationConvergenceMeasure(
                int dataId,
                bool suffices,
                int minimumIterationCount
                );

            virtual void measure(
                vector & oldValues,
                vector & newValues
                );

            virtual void newMeasurementSeries();

            virtual bool isConvergence();

            virtual void printState();

            const int minimumIterationCount;
            int currentIteration;
    };
}

#endif
