
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef GillebaartR3dFunction_H
#define GillebaartR3dFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class GillebaartR3dFunction : public RBFFunctionInterface
    {
public:

        GillebaartR3dFunction( double radius );

        virtual ~GillebaartR3dFunction();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual double correctionPowerCoeff();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
