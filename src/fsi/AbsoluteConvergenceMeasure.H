
/*
 * Copyright [2016] <David Blom>
 */

#ifndef AbsoluteConvergenceMeasure_H
#define AbsoluteConvergenceMeasure_H

#include "ConvergenceMeasure.H"
#include "fvCFD.H"

namespace fsi
{
    class AbsoluteConvergenceMeasure : public ConvergenceMeasure
    {
        public:
            AbsoluteConvergenceMeasure(
                int dataId,
                bool suffices,
                scalar convergenceLimit
                );

            virtual void measure(
                vector & oldValues,
                vector & newValues
                );

            virtual void newMeasurementSeries();

            virtual bool isConvergence();

            virtual void printState();

            scalar convergenceLimit;
            scalar normDiff;
    };
} // namespace fsi

#endif
