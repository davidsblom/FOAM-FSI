
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "GillebaartR3cFunction.H"

namespace rbf
{
    GillebaartR3cFunction::GillebaartR3cFunction( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    GillebaartR3cFunction::~GillebaartR3cFunction()
    {}

    double GillebaartR3cFunction::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        //return 1 - 10 * value * value * value + 15 * value * value * value * value - 6 * value * value * value * value * value;
        return 1 - 35 * value * value * value + 105 * value * value * value * value - 126 * value * value * value * value * value + 70 * value * value * value * value * value * value - 15 * value * value * value * value * value * value * value;
    }

    double GillebaartR3cFunction::correctionPower()
    {
        return 3;
    }

    double GillebaartR3cFunction::correctionPowerCoeff()
    {
        return 35;
    }


    void GillebaartR3cFunction::setRadius( double value )
    {
        assert( value > 0 );
        radius = value;
    }
}
