
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TPSFunction.H"
#include <iostream>

namespace rbf
{
    TPSFunction::TPSFunction()
    {}

    TPSFunction::~TPSFunction()
    {}

    double TPSFunction::evaluate( double value )
    {
        if ( value > 0 )
            return std::log( value ) * value * value;

        return 0;
    }

    double TPSFunction::correctionPower()
    {
        return 0;
    }

    void TPSFunction::setRadius( double value )
    {}
}
