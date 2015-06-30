
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "WendlandC4Function.H"

namespace rbf
{
    WendlandC4Function::WendlandC4Function( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    WendlandC4Function::~WendlandC4Function()
    {}

    double WendlandC4Function::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        //return std::pow( 1 - value, 6 ) * (35.0/3.0 * std::pow( value, 2 ) + 6.0 * value + 1.0);
        return ( 1 - value )*( 1 - value )*( 1 - value )*( 1 - value )*( 1 - value )*( 1 - value ) * (35.0/3.0 * value*value + 6.0 * value + 1);
    }
}
