
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC0Function_H
#define WendlandC0Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC0Function : public RBFFunctionInterface
    {
public:

        WendlandC0Function( double radius );

        virtual ~WendlandC0Function();

        virtual double evaluate( double value );

        virtual double correctionPower();

        virtual void setRadius( double value );

        double radius;
    };
}

#endif
