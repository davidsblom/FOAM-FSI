
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC2Function_H
#define WendlandC2Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC2Function : public RBFFunctionInterface
    {
public:

        WendlandC2Function( double radius );

        virtual ~WendlandC2Function();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual void setRadius( double value);

        double radius;
    };
}

#endif
