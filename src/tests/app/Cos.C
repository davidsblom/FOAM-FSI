
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "Cos.H"

using namespace sdc;

Cos::Cos(
    int nbTimeSteps,
    scalar dt,
    scalar endTime,
    scalar amplitude,
    scalar frequency
    )
    :
    f( 0 ),
    sol( 0 ),
    nbTimeSteps( nbTimeSteps ),
    dt( dt ),
    endTime( endTime ),
    amplitude( amplitude ),
    frequency( frequency ),
    t( 0 )
{
    assert( nbTimeSteps > 0 );
    assert( dt > 0 );
    assert( endTime > 0 );
    assert( amplitude > 0 );
    assert( frequency > 0 );

    sol = amplitude * ( 0.5 - 0.5 * std::cos( M_PI * frequency * t ) );
}

void Cos::evaluateFunction(
    const int,
    const fsi::vector &,
    const scalar t,
    fsi::vector & f
    )
{
    if ( t > 0 )
        f( 0 ) = this->f;
    else
        f( 0 ) = 0.5 * amplitude * std::sin( M_PI * frequency * t ) * M_PI * frequency;
}

void Cos::finalizeTimeStep()
{}

int Cos::getDOF()
{
    return 1;
}

void Cos::getSolution(
    fsi::vector & solution,
    fsi::vector &
    )
{
    solution( 0 ) = sol;
}

void Cos::setSolution(
    const fsi::vector &,
    const fsi::vector &
    )
{
    assert( false );
}

scalar Cos::getEndTime()
{
    return endTime;
}

scalar Cos::getTimeStep()
{
    return dt;
}

void Cos::nextTimeStep()
{}

void Cos::initTimeStep()
{}

void Cos::setNumberOfImplicitStages( int )
{}

void Cos::implicitSolve(
    bool,
    const int,
    const int,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs,
    fsi::vector & f,
    fsi::vector & result
    )
{
    f( 0 ) = 0.5 * amplitude * std::sin( M_PI * frequency * t ) * M_PI * frequency;
    result = dt * f + qold + rhs;

    result( 0 ) = amplitude * ( 0.5 - 0.5 * std::cos( M_PI * frequency * t ) );
    f = (result - qold - rhs) / dt;

    sol = result( 0 );

    this->f = f( 0 );
    this->t = t;
}
