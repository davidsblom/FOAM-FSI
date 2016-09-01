
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AdaptiveTimeStepper.H"
#include "PstreamReduceOps.H"

using namespace sdc;

AdaptiveTimeStepper::AdaptiveTimeStepper( bool enabled )
    :
    enabled( enabled ),
    filter( "h211b" ),
    tol( 1.0e-4 ),
    safetyFactor( 0.8 ),
    k( 0 ),
    cerrold( 0 ),
    rhoold( 0 ),
    timeStepIndex( 0 ),
    accepted( false ),
    previousTimeStepRejected( false ),
    endTime( 0 ),
    t( 0 )
{}

AdaptiveTimeStepper::AdaptiveTimeStepper(
    bool enabled,
    std::string filter,
    scalar tol,
    scalar safetyFactor
    )
    :
    enabled( enabled ),
    filter( filter ),
    tol( tol ),
    safetyFactor( safetyFactor ),
    k( 0 ),
    cerrold( 0 ),
    rhoold( 0 ),
    timeStepIndex( 0 ),
    accepted( false ),
    previousTimeStepRejected( false ),
    endTime( 0 ),
    t( 0 )
{
    assert( filter == "h211b" || filter == "elementary" || filter == "pi42" );
    assert( tol > 0 );
    assert( tol < 1 );
    assert( safetyFactor > 0 );
    assert( safetyFactor <= 1 );
}

AdaptiveTimeStepper::~AdaptiveTimeStepper()
{}

bool AdaptiveTimeStepper::determineNewTimeStep(
    const fsi::vector & errorEstimate,
    const fsi::vector & solution,
    const scalar computedTimeStep,
    scalar & newTimeStep
    )
{
    assert( endTime > 0 );

    scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );
    squaredNorm[Pstream::myProcNo()] = errorEstimate.squaredNorm();
    reduce( squaredNorm, sumOp<scalarList>() );
    scalar error = std::sqrt( sum( squaredNorm ) );

    squaredNorm = 0;
    squaredNorm[Pstream::myProcNo()] = solution.squaredNorm();
    reduce( squaredNorm, sumOp<scalarList>() );
    error /= std::sqrt( sum( squaredNorm ) );

    return determineNewTimeStep( error, computedTimeStep, newTimeStep );
}

bool AdaptiveTimeStepper::determineNewTimeStep(
    const scalar errorEstimate,
    const scalar computedTimeStep,
    scalar & newTimeStep
    )
{
    assert( endTime > 0 );
    assert( k > 0 );
    assert( errorEstimate > 0 );
    assert( enabled );

    scalar cerr = safetyFactor * tol / errorEstimate;
    accepted = false;
    scalar rho = 0;

    if ( tol / errorEstimate > 1 )
        accepted = true;

    if ( timeStepIndex == 0 || previousTimeStepRejected )
        rho = elementary( cerr );
    else
    {
        if ( filter == "h211b" )
            rho = h211b( cerr, cerrold, rhoold );

        if ( filter == "pi42" )
            rho = pi42( cerr, cerrold );
    }

    scalar ratio = limit( rho );

    newTimeStep = ratio * computedTimeStep;

    cerrold = cerr;
    rhoold = rho;

    if ( accepted )
    {
        previousTimeStepRejected = false;
        timeStepIndex++;
        t += computedTimeStep;

        if ( t + newTimeStep > endTime )
            newTimeStep = endTime - t;
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

scalar AdaptiveTimeStepper::elementary( const scalar c1 )
{
    // Elementary integrating controller for startup
    return std::pow( c1, 1.0 / k );
}

scalar AdaptiveTimeStepper::h211b(
    const scalar c1,
    const scalar c0,
    const scalar rho
    )
{
    // H211b digital filter, with standard parameter setting b = 4
    scalar b = 4;

    return std::pow( c1, 1.0 / b / k ) * std::pow( c0, 1.0 / b / k ) * std::pow( rho, -1.0 / b );
}

bool AdaptiveTimeStepper::isAccepted() const
{
    if ( not enabled )
        return true;

    return accepted;
}

bool AdaptiveTimeStepper::isEnabled() const
{
    return enabled;
}

bool AdaptiveTimeStepper::isPreviousStepAccepted() const
{
    return not previousTimeStepRejected;
}

scalar AdaptiveTimeStepper::limit( const scalar u )
{
    scalar kappa = 1;

    return 1.0 + kappa * std::atan( (u - 1) / kappa );
}

scalar AdaptiveTimeStepper::pi42(
    const scalar c1,
    const scalar c0
    )
{
    // PI.4.2 controller
    return std::pow( c1, 3.0 / 5.0 / k ) * std::pow( c0, -1.0 / 5.0 / k );
}

void AdaptiveTimeStepper::setEndTime( scalar endTime )
{
    assert( endTime > 0 );
    this->endTime = endTime;
}

void AdaptiveTimeStepper::setOrderEmbeddedMethod( int order )
{
    k = order + 1;
}
