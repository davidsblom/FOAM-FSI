
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearSolidSolver.H"
#include <boost/math/constants/constants.hpp>

using namespace tubeflow;

TubeFlowLinearSolidSolver::TubeFlowLinearSolidSolver(
    int N,
    scalar nu,
    scalar rho,
    scalar h,
    scalar L,
    scalar dt,
    scalar G,
    scalar E0,
    scalar r0,
    scalar T
    )
    :
    TubeFlowLinearizedSolidSolver( N, nu, rho, h, L, dt, G, E0, r0 ),
    T( T )
{
    assert( T > 0 );
    rn.fill( r0 );
    r.fill( r0 );
}

TubeFlowLinearSolidSolver::~TubeFlowLinearSolidSolver()
{}

void TubeFlowLinearSolidSolver::solve(
    const fsi::vector & p,
    fsi::vector & a
    )
{
    std::cout << "Solve solid domain" << std::endl;

    // Construct right hand size of linear system

    fsi::vector b( 2 * N ), x( 2 * N );

    for ( int i = 0; i < N; i++ )
    {
        b( i ) = rn( i );
        b( i + N ) = un( i ) + p( i ) * dt / (rho * h)
            + E0 * dt / (r0 * rho) / (1L - nu * nu);
    }

    b += rhs;

    // Solve for x

    x = lu.solve( b );

    // Retrieve solution

    r = x.head( N );
    u = x.tail( N );

    // Return area a
    a = boost::math::constants::pi<scalar>() * r.array() * r.array();

    data.col( 0 ) = a;
    this->p = p;
}
