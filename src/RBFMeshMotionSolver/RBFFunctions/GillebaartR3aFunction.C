
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "GillebaartR3aFunction.H"

namespace rbf
{
    GillebaartR3aFunction::GillebaartR3aFunction( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    GillebaartR3aFunction::~GillebaartR3aFunction()
    {}

    double GillebaartR3aFunction::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return 1 - 10 * value * value * value + 15 * value * value * value * value - 6 * value * value * value * value * value;
    }

    double GillebaartR3aFunction::correctionPower()
    {
        return 3;
    }

    double GillebaartR3aFunction::correctionPowerCoeff()
    {
        return 10;
    }

    void GillebaartR3aFunction::setRadius( double value )
    {
        assert( value > 0 );
        radius = value;
    }
}
