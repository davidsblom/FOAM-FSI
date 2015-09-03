
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "GillebaartR3Function.H"

namespace rbf
{
    GillebaartR3Function::GillebaartR3Function( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    GillebaartR3Function::~GillebaartR3Function()
    {}

    double GillebaartR3Function::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return 1 - 4 * value * value * value + 3 * value * value * value * value;
    }

    double GillebaartR3Function::correctionPower()
    {
        return 3;
    }

    double GillebaartR3Function::correctionPowerCoeff()
    {
        return 4;
    }

    void GillebaartR3Function::setRadius( double value )
    {
        assert( value > 0 );
        radius = value;
    }
}
