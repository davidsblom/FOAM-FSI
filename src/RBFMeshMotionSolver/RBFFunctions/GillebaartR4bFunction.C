
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "GillebaartR4bFunction.H"

namespace rbf
{
    GillebaartR4bFunction::GillebaartR4bFunction( double radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    GillebaartR4bFunction::~GillebaartR4bFunction()
    {}

    double GillebaartR4bFunction::evaluate( double value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        //return 1 - 10 * value * value * value + 15 * value * value * value * value - 6 * value * value * value * value * value;
        return 1 - 35 * value * value * value * value + 84 * value * value * value * value * value - 70 * value * value * value * value * value * value + 20 * value * value * value * value * value * value * value;
    }

    double GillebaartR4bFunction::correctionPower()
    {
        return 4;
    }

    double GillebaartR4bFunction::correctionPowerCoeff()
    {
        return 35;
    }


    void GillebaartR4bFunction::setRadius( double value )
    {
        assert( value > 0 );
        radius = value;
    }
}
