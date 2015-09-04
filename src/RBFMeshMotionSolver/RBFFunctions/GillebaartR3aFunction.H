
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef GillebaartR3aFunction_H
#define GillebaartR3aFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class GillebaartR3aFunction : public RBFFunctionInterface
    {
public:

        GillebaartR3aFunction( double radius );

        virtual ~GillebaartR3aFunction();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual double correctionPowerCoeff();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
