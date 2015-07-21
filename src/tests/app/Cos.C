
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "Cos.H"

using namespace sdc;

Cos::Cos(
    int nbTimeSteps,
    double dt,
    double endTime,
    double amplitude,
    double frequency
    )
    :
    f( 0 ),
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
    const int k,
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
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

double Cos::getScalingFactor()
{
    return 1.0;
}

void Cos::getSolution( Eigen::VectorXd & solution )
{
    solution( 0 ) = sol;
}

void Cos::setSolution(
    const Eigen::VectorXd & solution,
    const Eigen::VectorXd & f
    )
{
    assert( false );
}

double Cos::getEndTime()
{
    return endTime;
}

double Cos::getTimeStep()
{
    return dt;
}

void Cos::nextTimeStep()
{}

void Cos::initTimeStep()
{}

bool Cos::isRunning()
{
    if ( std::abs( t - endTime ) < 1.0e-13 )
        return false;

    return t <= endTime;
}

void Cos::setNumberOfImplicitStages( int k )
{}

void Cos::implicitSolve(
    bool corrector,
    const int k,
    const int kold,
    const double t,
    const double dt,
    const Eigen::VectorXd & qold,
    const Eigen::VectorXd & rhs,
    Eigen::VectorXd & f,
    Eigen::VectorXd & result
    )
{
    result( 0 ) = amplitude * ( 0.5 - 0.5 * std::cos( M_PI * frequency * t ) );

    f = (result - qold - rhs) / dt;

    sol = result( 0 );

    this->f = f( 0 );
}
