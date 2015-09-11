
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "GillebaartR3dFunction.H"

namespace rbf
{
    GillebaartR3dFunction::GillebaartR3dFunction( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    GillebaartR3dFunction::~GillebaartR3dFunction()
    {}

    double GillebaartR3dFunction::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return 1 - 56 * value * value * value + 210 * value * value * value * value - 336 * value * value * value * value * value + 280 * value * value * value * value * value * value - 120 * value * value * value * value * value * value * value + 21 * value * value * value * value * value * value * value * value;
    }

    double GillebaartR3dFunction::correctionPower()
    {
        return 3;
    }

    double GillebaartR3dFunction::correctionPowerCoeff()
    {
        return 56;
    }


    void GillebaartR3dFunction::setRadius( double value )
    {
        assert( value > 0 );
        radius = value;
    }
}
