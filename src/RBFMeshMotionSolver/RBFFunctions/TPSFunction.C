
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TPSFunction.H"

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
}
