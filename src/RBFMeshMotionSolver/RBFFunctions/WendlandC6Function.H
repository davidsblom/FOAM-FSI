
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC6Function_H
#define WendlandC6Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC6Function : public RBFFunctionInterface
    {
public:

        WendlandC6Function( double radius );

        virtual ~WendlandC6Function();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual double correctionPowerCoeff();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
