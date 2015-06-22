
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AdaptiveTimeStepper.H"

using namespace sdc;

AdaptiveTimeStepper::AdaptiveTimeStepper(
    std::string filter,
    double tol,
    int k
    )
    :
    filter( filter ),
    tol( tol ),
    k( k ),
    cerrold( 0 ),
    rhoold( 0 ),
    timeStepIndex( 0 )
{
    assert( filter == "h211b" || filter == "elementary" || filter == "pi42" );
    assert( tol > 0 );
    assert( tol < 1 );
    assert( k > 0 );
}

AdaptiveTimeStepper::~AdaptiveTimeStepper()
{}

bool AdaptiveTimeStepper::determineNewTimeStep(
    const double errorEstimate,
    const double computedTimeStep,
    double & newTimeStep
    )
{
    double cerr = tol / errorEstimate;
    bool accepted = false;
    double rho = 0;

    if ( cerr > 1 )
        accepted = true;

    if ( timeStepIndex == 0 )
        rho = elementary( cerr, cerrold, rhoold );

    if ( timeStepIndex > 0 )
    {
        if ( filter == "h211b" )
            rho = h211b( cerr, cerrold, rhoold );

        if ( filter == "pi42" )
            rho = pi42( cerr, cerrold, rhoold );
    }

    double ratio = limit( rho );

    newTimeStep = ratio * computedTimeStep;

    cerrold = cerr;
    rhoold = rho;

    if ( accepted )
        timeStepIndex++;

    return accepted;
}

double AdaptiveTimeStepper::elementary(
    const double c1,
    const double c0,
    const double rho
    )
{
    // Elementary integrating controller for startup
    return pow( c1, 1.0 / k );
}

double AdaptiveTimeStepper::h211b(
    const double c1,
    const double c0,
    const double rho
    )
{
    // H211b digital filter, with standard parameter setting b = 4
    double b = 4;

    return pow( c1, 1.0 / b / k ) * pow( c0, 1.0 / b / k ) * pow( rho, -1.0 / b );
}

double AdaptiveTimeStepper::limit( const double u )
{
    double kappa = 1;

    return 1.0 + kappa * atan( (u - 1) / kappa );
}

double AdaptiveTimeStepper::pi42(
    const double c1,
    const double c0,
    const double rho
    )
{
    // PI.4.2 controller
    return pow( c1, 3.0 / 5.0 / k ) * pow( c0, -1.0 / 5.0 / k );
}
