
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "WendlandC0Function.H"

namespace rbf
{
    WendlandC0Function::WendlandC0Function( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    WendlandC0Function::~WendlandC0Function()
    {}

    double WendlandC0Function::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return std::pow( 1 - value, 2 );
    }
}
