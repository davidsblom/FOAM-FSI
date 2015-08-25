
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ConvergenceMeasure_H
#define ConvergenceMeasure_H

#include "BaseMultiLevelSolver.H"
#include <iomanip>

namespace fsi
{
    class ConvergenceMeasure
    {
public:

        ConvergenceMeasure( int dataId )
            :
            dataId( dataId ),
            isConvergence_( false )
        {}

        virtual ~ConvergenceMeasure(){}

        virtual void measure(
            vector & oldValues,
            vector & newValues
            ) = 0;

        virtual void newMeasurementSeries() = 0;

        virtual bool isConvergence() = 0;

        virtual void printState() = 0;

        const int dataId;

        bool isConvergence_;
    };
}

#endif
