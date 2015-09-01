
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC4Function_H
#define WendlandC4Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC4Function : public RBFFunctionInterface
    {
public:

        WendlandC4Function( double radius );

        virtual ~WendlandC4Function();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
