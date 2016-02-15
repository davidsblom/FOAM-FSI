
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowLinearizedFluidSolver.H"

using namespace tubeflow;

TubeFlowLinearizedFluidSolver::TubeFlowLinearizedFluidSolver(
    int N,
    scalar p0,
    scalar r0,
    scalar u0,
    scalar rho,
    scalar E,
    scalar h,
    scalar T,
    scalar dt,
    scalar L
    ) :
    BaseMultiLevelSolver( N, 1, p0 ),
    p0( p0 ),
    r0( r0 ),
    u0( u0 ),
    a0( M_PI * r0 * r0 ),
    rho( rho ),
    h( h ),
    T( T ),
    dt( dt ),
    dx( L / N ),
    tau( u0 * dt / L ),
    cmk( std::sqrt( E * h / (2 * rho * r0) ) ),
    alpha( a0 / (u0 + dx / dt) ),
    grid(),
    un( N ),
    an( N ),
    pn( N ),
    u( N ),
    a( N ),
    p( N ),
    lu()
{
    assert( N > 0 );
    assert( p0 >= 0 );
    assert( r0 > 0 );
    assert( T > 0 );
    assert( T >= dt );
    assert( rho > 0 );
    assert( u0 >= 0 );
    assert( cmk > 0 );
    assert( h > 0 );

    p.setZero();
    pn.setZero();
    u.setZero();
    un.setZero();
    a.setZero();
    an.setZero();
    data.fill( p0 );

    factorizeMatrix();
}

TubeFlowLinearizedFluidSolver::~TubeFlowLinearizedFluidSolver()
{}

void TubeFlowLinearizedFluidSolver::calcGrid()
{
    if ( grid.rows() != N )
    {
        grid.resize( N, 1 );
        grid.setZero();

        for ( int i = 0; i < N; i++ )
            grid( i, 0 ) = dx * i + 0.5 * dx;
    }
}

void TubeFlowLinearizedFluidSolver::factorizeMatrix()
{
    // Initialize matrix A

    matrix A( 2 * N, 2 * N );
    A.setZero();

    // Build A matrix

    for ( int i = 0; i < N; i++ )
    {
        if ( i == 0 )
        {
            // Inlet

            // Velocity inlet boundary condition
            A( i, i ) = 1;

            // Pressure inlet boundary condition
            A( i + N, i + N ) = 1;
            A( i + N, i + N + 1 ) = -2;
            A( i + N, i + N + 2 ) = 1;
        }

        if ( i == N - 1 )
        {
            // Outlet

            // Velocity outlet boundary condition
            A( i, i ) = 1;
            A( i, i - 1 ) = -2;
            A( i, i - 2 ) = 1;

            // Pressure outlet boundary condition
            A( i + N, i + N ) = 1;
            A( i + N, i ) = -cmk * rho;
        }

        if ( i > 0 && i < N - 1 )
        {
            // Continuity equation

            A( i, i + 1 ) = 0.5 * a0;
            A( i, i - 1 ) = -0.5 * a0;
            A( i, i + N ) = 2 * alpha / rho;
            A( i, i + N + 1 ) = -alpha / rho;
            A( i, i + N - 1 ) = -alpha / rho;

            // Momentum equation

            A( i + N, i ) = a0 * dx / dt + u0 * a0;
            A( i + N, i + 1 ) = 0.5 * u0 * a0;
            A( i + N, i - 1 ) = -1.5 * u0 * a0;
            A( i + N, i + N + 1 ) = a0 / (2 * rho);
            A( i + N, i + N - 1 ) = -a0 / (2 * rho);
        }
    }

    // Compute LU decomposition

    lu.compute( A );
}

void TubeFlowLinearizedFluidSolver::finalizeTimeStep()
{
    assert( init );

    un = u;
    pn = p;
    an = a;

    init = false;
}

void TubeFlowLinearizedFluidSolver::getReadPositions( matrix & readPositions )
{
    calcGrid();

    readPositions = grid;
}

void TubeFlowLinearizedFluidSolver::getWritePositions( matrix & writePositions )
{
    calcGrid();

    writePositions = grid;
}

void TubeFlowLinearizedFluidSolver::initTimeStep()
{
    assert( !init );

    timeIndex++;
    t = timeIndex * dt;

    init = true;
}

bool TubeFlowLinearizedFluidSolver::isRunning()
{
    assert( !init );

    return t < T;
}

void TubeFlowLinearizedFluidSolver::resetSolution()
{}

void TubeFlowLinearizedFluidSolver::solve(
    const matrix & input,
    matrix & output
    )
{
    assert( input.rows() == N );
    assert( input.cols() == 1 );

    // Map the matrices to fsi::vectors
    fsi::vector a = Eigen::Map<const fsi::vector> ( input.data(), input.rows() );
    fsi::vector p;

    // Call solve with the correct arguments
    solve( a, p );

    output = p;

    assert( output.rows() == N );
    assert( output.cols() == 1 );
}

void TubeFlowLinearizedFluidSolver::solve(
    const fsi::vector & a,
    fsi::vector & p
    )
{
    this->a = a.array() - a0;

    std::cout << "Solve fluid domain" << std::endl;

    // Construct right hand size of linear system

    fsi::vector b( 2 * N ), x( 2 * N );
    b.setZero();

    for ( int i = 1; i < N - 1; i++ )
    {
        // Continuity equation rhs

        b( i ) = dx / dt * ( an( i ) - this->a( i ) );
        b( i ) += 0.5 * u0 * ( this->a( i - 1 ) - this->a( i + 1 ) );

        // Momentum equation rhs

        b( i + N ) = u0 * dx / dt * ( an( i ) - this->a( i ) );
        b( i + N ) += 0.5 * u0 * u0 * ( this->a( i - 1 ) - this->a( i + 1 ) );
        b( i + N ) += a0 * dx / dt * un( i );
    }

    // Boundary conditions

    // Velocity inlet condition
    b( 0 ) = u0 / 10.0 * std::pow( std::sin( M_PI * timeIndex * tau ), 2 );

    // Pressure inlet
    b( N ) = 0;

    // Velocity outlet condition
    b( N - 1 ) = 0;

    // Pressure outlet condition
    b( 2 * N - 1 ) = -cmk * rho * un( N - 1 ) + pn( N - 1 );

    // Solve for x

    x = lu.solve( b );

    // Retrieve solution

    u = x.head( N );
    this->p = x.tail( N );

    // Return pressure p
    p = this->p.array() + p0;

    data.col( 0 ) = p;
}
