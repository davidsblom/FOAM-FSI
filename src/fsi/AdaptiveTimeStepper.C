
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AdaptiveTimeStepper.H"

using namespace sdc;

AdaptiveTimeStepper::AdaptiveTimeStepper(
    bool enabled,
    std::string filter,
    double tol,
    double safetyFactor,
    int k
    )
    :
    enabled( enabled ),
    filter( filter ),
    tol( tol ),
    safetyFactor( safetyFactor ),
    k( k ),
    cerrold( 0 ),
    rhoold( 0 ),
    timeStepIndex( 0 ),
    accepted( false ),
    previousTimeStepRejected( false )
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
    if ( not enabled )
    {
        newTimeStep = computedTimeStep;
        return true;
    }

    double cerr = safetyFactor * tol / errorEstimate;
    accepted = false;
    double rho = 0;

    if ( tol / errorEstimate > 1 )
        accepted = true;

    if ( timeStepIndex == 0 || previousTimeStepRejected )
        rho = elementary( cerr, cerrold, rhoold );
    else
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
    {
        previousTimeStepRejected = false;
        timeStepIndex++;
    }
    else
        previousTimeStepRejected = true;

    Info << "adaptive time: error = " << errorEstimate;
    Info << ", accepted = ";

    if ( accepted )
        Info << "true";
    else
        Info << "false";

    Info << ", current time step = " << computedTimeStep;
    Info << ", new time step = " << newTimeStep << endl;

    return accepted;
}

double AdaptiveTimeStepper::elementary(
    const double c1,
    const double c0,
    const double rho
    )
{
    // Elementary integrating controller for startup
    return std::pow( c1, 1.0 / k );
}

double AdaptiveTimeStepper::h211b(
    const double c1,
    const double c0,
    const double rho
    )
{
    // H211b digital filter, with standard parameter setting b = 4
    double b = 4;

    return std::pow( c1, 1.0 / b / k ) * std::pow( c0, 1.0 / b / k ) * std::pow( rho, -1.0 / b );
}

bool AdaptiveTimeStepper::isAccepted()
{
    if ( not enabled )
        return true;

    return accepted;
}

bool AdaptiveTimeStepper::isEnabled()
{
    return enabled;
}

bool AdaptiveTimeStepper::isPreviousStepAccepted()
{
    return not previousTimeStepRejected;
}

double AdaptiveTimeStepper::limit( const double u )
{
    double kappa = 1;

    return 1.0 + kappa * std::atan( (u - 1) / kappa );
}

double AdaptiveTimeStepper::pi42(
    const double c1,
    const double c0,
    const double rho
    )
{
    // PI.4.2 controller
    return std::pow( c1, 3.0 / 5.0 / k ) * std::pow( c0, -1.0 / 5.0 / k );
}
