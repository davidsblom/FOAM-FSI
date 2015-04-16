
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "WendlandC2Function.H"

namespace rbf
{
    WendlandC2Function::WendlandC2Function( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    WendlandC2Function::~WendlandC2Function()
    {}

    double WendlandC2Function::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return std::pow( 1 - value, 4 ) * (4 * value + 1);
    }
}
