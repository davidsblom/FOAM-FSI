
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef GillebaartR4bFunction_H
#define GillebaartR4bFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class GillebaartR4bFunction : public RBFFunctionInterface
    {
public:

        GillebaartR4bFunction( double radius );

        virtual ~GillebaartR4bFunction();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual double correctionPowerCoeff();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
