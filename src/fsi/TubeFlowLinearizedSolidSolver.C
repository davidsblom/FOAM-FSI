
/*
 * Copyright [2016] <David Blom>
 */

#include "TubeFlowLinearizedSolidSolver.H"
#include <boost/math/constants/constants.hpp>

namespace tubeflow {

TubeFlowLinearizedSolidSolver::TubeFlowLinearizedSolidSolver(
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
    fsi::BaseMultiLevelSolver( N, 1, boost::math::constants::pi<scalar>() * r0 * r0 ),
    N( N ),
    timeOrder( 1 ),
    dt( dt ),
    r0( r0 ),
    kappa( 2.0L * (1.0L + nu) / (4.0L + 3.0L * nu) ),
    dx( L / N ),
    G( G ),
    E0( E0 ),
    nu( nu ),
    h( h ),
    rho( rho ),
    T( T ),
    grid(),
    un( N ),
    rn( N ),
    u( N ),
    r( N ),
    rhs( 2 * N ),
    p( N ),
    lu( nullptr ),
    alpha(
{
    -1., 1.
} ),
    beta( 1. ),
    rStages(),
    uStages()
{
    assert( N > 0 );
    assert( nu > 0 );
    assert( nu < 1 );
    assert( h > 0 );
    assert( kappa > 0 );
    assert( L > 0 );
    assert( G > 0 );
    assert( r0 > 0 );
    assert( E0 > 0 );
    assert( T > 0 );

    un.setZero();
    u.setZero();
    rn.setZero();
    r.setZero();
    rhs.setZero();
    p.setZero();

    scalar a0 = boost::math::constants::pi<scalar>() * r0 * r0;
    data.fill( a0 );

    factorizeMatrix();

    rStages.push_back( rn );
    rStages.push_back( r );
    uStages.push_back( un );
    uStages.push_back( u );
}

TubeFlowLinearizedSolidSolver::TubeFlowLinearizedSolidSolver(
    int N,
    scalar nu,
    scalar rho,
    scalar h,
    scalar L,
    scalar dt,
    scalar G,
    scalar E0,
    scalar r0,
    scalar T,
    int timeOrder
    )
    :
    BaseMultiLevelSolver( N, 1, boost::math::constants::pi<scalar>() * r0 * r0 ),
    N( N ),
    timeOrder( timeOrder ),
    dt( dt ),
    r0( r0 ),
    kappa( 2.0L * (1.0L + nu) / (4.0L + 3.0L * nu) ),
    dx( L / N ),
    G( G ),
    E0( E0 ),
    nu( nu ),
    h( h ),
    rho( rho ),
    T( T ),
    grid(),
    un( N ),
    rn( N ),
    u( N ),
    r( N ),
    rhs( 2 * N ),
    p( N ),
    lu( nullptr ),
    alpha(
{
    -1, 1
} ),
    beta( 1. )
{
    assert( N > 0 );
    assert( nu > 0 );
    assert( nu < 1 );
    assert( h > 0 );
    assert( kappa > 0 );
    assert( L > 0 );
    assert( G > 0 );
    assert( r0 > 0 );
    assert( E0 > 0 );
    assert( T > 0 );
    assert( timeOrder == 1 || timeOrder == 2 );

    un.setZero();
    u.setZero();
    rn.setZero();
    r.setZero();
    rhs.setZero();
    p.setZero();

    scalar a0 = boost::math::constants::pi<scalar>() * r0 * r0;
    data.fill( a0 );

    factorizeMatrix();

    rStages.push_back( rn );
    rStages.push_back( r );
    uStages.push_back( un );
    uStages.push_back( u );
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
            grid( i, 0 ) = dx * i + 0.5L * dx;
    }
}

void TubeFlowLinearizedSolidSolver::factorizeMatrix()
{
    // Construct matrix A

    fsi::matrix A( 2 * N, 2 * N );
    A.setZero();

    // Build A matrix

    for ( int i = 0; i < N; i++ )
    {
        A( i, i ) = alpha.back();
        A( i, i + N ) = -dt * beta;

        if ( i == 0 )
        {
            // Boundary conditions at inlet
            // r(0) = r(1)

            A( i + N, i ) = kappa * G * dt * beta / (rho * dx * dx) + E0 * dt * beta / (1L - nu * nu) / (rho * r0 * r0);
            A( i + N, i + 1 ) = -kappa * G * dt * beta / (rho * dx * dx);

            A( i + N, i + N ) = alpha.back();
        }

        if ( i == N - 1 )
        {
            // Boundary conditions at outlet
            // r(N) = r(N-1)

            A( i + N, i ) = kappa * G * dt * beta / (rho * dx * dx) + E0 * dt * beta / (1L - nu * nu) / (rho * r0 * r0);
            A( i + N, i - 1 ) = -kappa * G * dt * beta / (rho * dx * dx);

            A( i + N, i + N ) = alpha.back();
        }

        if ( i > 0 && i < N - 1 )
        {
            A( i + N, i ) = 2L * kappa * G * dt * beta / (rho * dx * dx) + E0 * dt * beta / (1L - nu * nu) / (rho * r0 * r0);
            A( i + N, i + 1 ) = -kappa * G * dt * beta / (rho * dx * dx);
            A( i + N, i - 1 ) = -kappa * G * dt * beta / (rho * dx * dx);

            A( i + N, i + N ) = alpha.back();
        }
    }

    lu = std::shared_ptr<Eigen::FullPivLU<fsi::matrix> > ( new Eigen::FullPivLU<fsi::matrix>( A ) );
}

void TubeFlowLinearizedSolidSolver::finalizeTimeStep()
{
    assert( init );

    un = u;
    rn = r;

    rStages.push_back( r );
    uStages.push_back( u );

    if ( timeOrder == 2 && alpha.size() == 2 )
    {
        alpha = {
            1. / 3., -4. / 3., 1.
        };
        beta = 2. / 3.;
        factorizeMatrix();
    }

    while ( rStages.size() > alpha.size() )
    {
        rStages.erase( rStages.begin() );
        uStages.erase( uStages.begin() );
    }

    assert( rStages.size() == uStages.size() );

    init = false;
}

void TubeFlowLinearizedSolidSolver::getReadPositions( fsi::matrix & readPositions )
{
    calcGrid();

    readPositions = grid;
}

void TubeFlowLinearizedSolidSolver::getWritePositions( fsi::matrix & writePositions )
{
    calcGrid();

    writePositions = grid;
}

void TubeFlowLinearizedSolidSolver::initTimeStep()
{
    assert( !init );

    timeIndex++;
    t = timeIndex * dt;

    init = true;
}

bool TubeFlowLinearizedSolidSolver::isRunning()
{
    assert( !init );

    return t < T;
}

void TubeFlowLinearizedSolidSolver::resetSolution()
{}

void TubeFlowLinearizedSolidSolver::solve(
    const fsi::matrix & input,
    fsi::matrix & output
    )
{
    assert( input.rows() == N );
    assert( input.cols() == 1 );
    assert( init );

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
    assert( init );
    assert( rStages.size() == uStages.size() );
    assert( rStages.size() == alpha.size() );

    std::cout << "Solve solid domain" << std::endl;

    // Construct right hand size of linear system

    fsi::vector b( 2 * N ), x( 2 * N );

    for ( int i = 0; i < N; i++ )
    {
        b( i ) = 0;

        for ( unsigned int k = 0; k < alpha.size() - 1; k++ )
            b( i ) -= alpha[k] * rStages[k]( i );

        b( i + N ) = p( i ) * dt * beta / (rho * h);

        for ( unsigned int k = 0; k < alpha.size() - 1; k++ )
            b( i + N ) -= alpha[k] * uStages[k]( i );
    }

    if ( rhs.rows() > 0 )
        b += rhs;

    // Solve for x

    x = lu->solve( b );

    // Retrieve solution

    r = x.head( N );
    u = x.tail( N );

    // Return area a
    a = r.array() + r0;
    a = boost::math::constants::pi<scalar>() * a.array() * a.array();

    data.col( 0 ) = a;
    this->p = p;
}

void TubeFlowLinearizedSolidSolver::run()
{
    while ( isRunning() )
        solveTimeStep();
}

void TubeFlowLinearizedSolidSolver::solveTimeStep()
{
    initTimeStep();
    fsi::vector a;
    solve( p, a );
    finalizeTimeStep();
}

} // namespace tubeflow
