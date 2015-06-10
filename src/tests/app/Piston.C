
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "Piston.H"

Piston::Piston(
    int nbTimeSteps,
    double dt,
    double q0,
    double qdot0,
    double As,
    double Ac,
    double omega
    )
    :
    SDCSolver(),
    nbTimeSteps( nbTimeSteps ),
    dt( dt ),
    q0( q0 ),
    c1( q0 + Ac / std::pow( omega, 2 ) ),
    c2( qdot0 + As / omega ),
    As( As ),
    Ac( Ac ),
    omega( omega ),
    N( 2 ),
    q( q0 ),
    qdot( qdot0 ),
    t( 0 ),
    timeIndex( 0 )
{
    assert( nbTimeSteps > 0 );
    assert( dt > 0 );
}

Piston::~Piston()
{}

double Piston::referenceSolution( double t )
{
    double result = c1;
    result += c2 * t;
    result += -As / std::pow( omega, 2 ) * std::sin( omega * t );
    result += -Ac / std::pow( omega, 2 ) * std::cos( omega * t );

    return result;
}

void Piston::evaluateFunction(
    const int k,
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
    )
{
    assert( f.rows() == 2 );

    f( 0 ) = As * std::sin( omega * t );
    f( 0 ) += Ac * std::cos( omega * t );
    f( 1 ) = q( 0 );
}

void Piston::nextTimeStep()
{
    timeIndex++;
}

int Piston::getDOF()
{
    return 2;
}

void Piston::getSolution( Eigen::VectorXd & solution )
{
    assert( solution.rows() == 2 );
    solution( 0 ) = qdot;
    solution( 1 ) = q;
}

double Piston::getTimeStep()
{
    return dt;
}

bool Piston::isRunning()
{
    return timeIndex < nbTimeSteps;
}

void Piston::run()
{
    Eigen::VectorXd q( nbTimeSteps + 1 ), qdot( nbTimeSteps + 1 ), qold( 2 ), f( 2 ), rhs( 2 ), result( 2 );

    q( 0 ) = -Ac;
    qdot( 0 ) = q( 0 );

    for ( int i = 1; i < nbTimeSteps + 1; i++ )
    {
        double t = dt * i;

        qold << qdot( i - 1 ), q( i - 1 );

        f.setZero();
        rhs.setZero();

        implicitSolve( false, 0, t, dt, 1, qold, rhs, f, result );

        qdot( i ) = result( 0 );
        q( i ) = result( 1 );
    }
}

void Piston::implicitSolve(
    bool corrector,
    const int k,
    const double t,
    const double dt,
    const double Akk,
    const Eigen::VectorXd & qold,
    const Eigen::VectorXd & rhs,
    Eigen::VectorXd & f,
    Eigen::VectorXd & result
    )
{
    assert( f.rows() == 2 );
    assert( rhs.rows() == 2 );
    assert( result.rows() == 2 );

    double deltaT = dt * Akk;

    f( 0 ) = As * std::sin( omega * t );
    f( 0 ) += Ac * std::cos( omega * t );

    // qdot
    result( 0 ) = qold( 0 ) + deltaT * f( 0 ) + rhs( 0 );

    // q
    result( 1 ) = qold( 1 ) + std::pow( deltaT, 2 ) * f( 0 ) + deltaT * qold( 0 ) + deltaT * rhs( 0 ) + rhs( 1 );

    f( 1 ) = result( 0 );

    qdot = result( 0 );
    q = result( 1 );

    this->t = t;
}
