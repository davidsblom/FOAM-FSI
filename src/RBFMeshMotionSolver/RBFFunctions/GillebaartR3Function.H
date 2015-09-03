
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef GillebaartR3Function_H
#define GillebaartR3Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class GillebaartR3Function : public RBFFunctionInterface
    {
public:

        GillebaartR3Function( double radius );

        virtual ~GillebaartR3Function();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual double correctionPowerCoeff();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
