
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

        return (1 - value) * (1 - value) * (1 - value) * (1 - value) * (4 * value + 1);
    }

    double WendlandC2Function::correctionPower()
    {
        return 2;
    }

    double WendlandC2Function::correctionPowerCoeff()
    {
        return 10;
    }


    void WendlandC2Function::setRadius( double value )
    {
        assert( value > 0 );
        radius = value;
    }
}
