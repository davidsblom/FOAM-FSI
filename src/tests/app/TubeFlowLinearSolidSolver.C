
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearSolidSolver.H"

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
    scalar r0
    )
    :
    TubeFlowLinearizedSolidSolver( N, nu, rho, h, L, dt, G, E0, r0 ),
    h( h ),
    E0( E0 ),
    nu( nu )
{
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
    b.setZero();

    for ( int i = 1; i < N - 1; i++ )
    {
        // Velocity equation based on backward Euler time stepping

        b( i ) = -rn( i );

        // Newton's equation rhs

        b( i + N ) = p( i ) + alpha * un( i ) + E0 * h / (1.0 - nu * nu) * 1.0 / r0;
    }

    // Boundary conditions

    b( 0 ) = 0;
    b( N - 1 ) = 0;
    b( N ) = 0;
    b( 2 * N - 1 ) = 0;

    // Solve for x

    x = lu.solve( b );

    // Retrieve solution

    u = x.head( N );
    r = x.tail( N );

    // Return area a
    a = M_PI * r.array() * r.array();

    data.col( 0 ) = a;
}
