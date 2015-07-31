
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearizedSolidSolver.H"

using namespace tubeflow;

TubeFlowLinearizedSolidSolver::TubeFlowLinearizedSolidSolver(
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
    BaseMultiLevelSolver( N, 1, M_PI * r0 * r0 ),
    N( N ),
    dt( dt ),
    r0( r0 ),
    kappa( 2.0 * (1.0 + nu) / (4.0 + 3.0 * nu) ),
    dx( L / N ),
    alpha( rho * h / dt ),
    beta( -kappa * G * h / (dx * dx) ),
    gamma( E0 * h / (1.0 - nu * nu) * ( 1.0 / (r0 * r0) ) ),
    grid(),
    un( N ),
    rn( N ),
    u( N ),
    r( N ),
    lu()
{
    assert( N > 0 );
    assert( nu > 0 );
    assert( nu < 1 );
    assert( h > 0 );
    assert( kappa > 0 );
    assert( alpha > 0 );
    assert( L > 0 );
    assert( G > 0 );
    assert( r0 > 0 );

    un.setZero();
    u.setZero();
    rn.setZero();
    r.setZero();

    scalar a0 = M_PI * r0 * r0;
    data.fill( a0 );

    factorizeMatrix();
}

TubeFlowLinearizedSolidSolver::~TubeFlowLinearizedSolidSolver()
{}

void TubeFlowLinearizedSolidSolver::calcGrid()
{
    if ( grid.rows() != N )
    {
        grid.resize( N, 1 );
        grid.setZero();

        for ( int i = 0; i < N; i++ )
            grid( i, 0 ) = dx * i + 0.5 * dx;
    }
}

void TubeFlowLinearizedSolidSolver::factorizeMatrix()
{
    // Construct matrix A

    matrix A( 2 * N, 2 * N );
    A.setZero();

    // Build A matrix

    for ( int i = 0; i < N; i++ )
    {
        // Velocity equation based on backward Euler time stepping

        // Newton's equation

        if ( i == 0 )
        {
            // Boundary conditions at inlet
            // r(0) = r(1)

            A( i, i ) = 1;

            A( i + N, i + N ) = 1;
            A( i + N, i + N + 1 ) = -1;
        }

        if ( i == N - 1 )
        {
            // Boundary conditions at outlet
            // r(N) = r(N-1)

            A( i, i ) = 1;

            A( i + N, i + N ) = 1;
            A( i + N, i + N - 1 ) = -1;
        }

        if ( i > 0 && i < N - 1 )
        {
            A( i, i ) = dt;
            A( i, i + N ) = -1;

            A( i + N, i ) = alpha;
            A( i + N, i + N + 1 ) = beta;
            A( i + N, i + N ) = gamma - 2 * beta;
            A( i + N, i + N - 1 ) = beta;
        }
    }

    lu.compute( A );
}

void TubeFlowLinearizedSolidSolver::finalizeTimeStep()
{
    assert( init );

    un = u;
    rn = r;

    init = false;
}

void TubeFlowLinearizedSolidSolver::getReadPositions( matrix & readPositions )
{
    calcGrid();

    readPositions = grid;
}

void TubeFlowLinearizedSolidSolver::getWritePositions( matrix & writePositions )
{
    calcGrid();

    writePositions = grid;
}

void TubeFlowLinearizedSolidSolver::initTimeStep()
{
    assert( !init );

    init = true;
}

bool TubeFlowLinearizedSolidSolver::isRunning()
{
    return true;
}

void TubeFlowLinearizedSolidSolver::resetSolution()
{}

void TubeFlowLinearizedSolidSolver::solve(
    const matrix & input,
    matrix & output
    )
{
    assert( input.rows() == N );
    assert( input.cols() == 1 );

    // Map the matrices to fsi::vectors
    fsi::vector p = Eigen::Map< const fsi::vector> ( input.data(), input.rows() );
    fsi::vector a;

    // Call solve with the correct arguments
    solve( p, a );

    output = a;

    assert( output.rows() == N );
    assert( output.cols() == 1 );
}

void TubeFlowLinearizedSolidSolver::solve(
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

        b( i + N ) = p( i ) + alpha * un( i );
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
    a = r.array() + r0;
    a = M_PI * a.array() * a.array();

    data.col( 0 ) = a;
}
