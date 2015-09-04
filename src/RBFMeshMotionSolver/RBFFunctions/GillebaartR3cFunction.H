
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef GillebaartR3cFunction_H
#define GillebaartR3cFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class GillebaartR3cFunction : public RBFFunctionInterface
    {
public:

        GillebaartR3cFunction( double radius );

        virtual ~GillebaartR3cFunction();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual double correctionPowerCoeff();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
