
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowExplicitLinearSolidSolver.H"
#include <boost/math/constants/constants.hpp>

using namespace tubeflow;

TubeFlowExplicitLinearSolidSolver::TubeFlowExplicitLinearSolidSolver(
    int N,
    scalar nu,
    scalar rho,
    scalar h,
    scalar L,
    scalar dt,
    scalar G,
    scalar E0,
    scalar r0,
    scalar p0,
    scalar T
    )
    :
    BaseMultiLevelSolver( N, 1, boost::math::constants::pi<scalar>() * r0 * r0 ),
    N( N ),
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
    pn( N )
{
    assert( N > 0 );
    assert( nu > 0 );
    assert( nu < 1 );
    assert( h > 0 );
    assert( kappa > 0 );
    assert( L > 0 );
    assert( G > 0 );
    assert( r0 > 0 );
    assert( T > 0 );

    un.setZero();
    u.setZero();
    rhs.setZero();
    p.fill( p0 );
    pn.fill( p0 );
    r.setZero();
    rn.setZero();

    scalar a0 = boost::math::constants::pi<scalar>() * r0 * r0;
    data.fill( a0 );
}

TubeFlowExplicitLinearSolidSolver::~TubeFlowExplicitLinearSolidSolver()
{}

void TubeFlowExplicitLinearSolidSolver::calcGrid()
{
    if ( grid.rows() != N )
    {
        grid.resize( N, 1 );
        grid.setZero();

        for ( int i = 0; i < N; i++ )
            grid( i, 0 ) = dx * i + 0.5L * dx;
    }
}

void TubeFlowExplicitLinearSolidSolver::finalizeTimeStep()
{
    assert( init );

    un = u;
    rn = r;
    pn = p;

    init = false;
}

void TubeFlowExplicitLinearSolidSolver::getReadPositions( matrix & readPositions )
{
    calcGrid();

    readPositions = grid;
}

void TubeFlowExplicitLinearSolidSolver::getWritePositions( matrix & writePositions )
{
    calcGrid();

    writePositions = grid;
}

void TubeFlowExplicitLinearSolidSolver::initTimeStep()
{
    assert( !init );

    init = true;
}

bool TubeFlowExplicitLinearSolidSolver::isRunning()
{
    return true;
}

void TubeFlowExplicitLinearSolidSolver::resetSolution()
{}

void TubeFlowExplicitLinearSolidSolver::solve(
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

void TubeFlowExplicitLinearSolidSolver::solve(
    const fsi::vector & p,
    fsi::vector & a
    )
{
    std::cout << "Solve solid domain" << std::endl;

    int nbStages = 4;
    std::vector<fsi::vector> kStages( nbStages );

    fsi::matrix A( 4, 4 );
    A.setZero();
    A( 1, 0 ) = 1. / 3.;
    A( 2, 0 ) = -1. / 3.;
    A( 2, 1 ) = 1.;
    A( 3, 0 ) = 1.;
    A( 3, 1 ) = -1.;
    A( 3, 2 ) = 1.;

    std::vector<scalar> b_RK4 = {
        1. / 8., 3. / 8., 3. / 8., 1. / 8.
    };

    for ( int iStage = 0; iStage < nbStages; iStage++ )
    {
        fsi::vector rStage = rn, uStage = un;

        for ( int i = 0; i < iStage; i++ )
        {
            rStage.array() += dt * A( iStage, i ) * kStages[i].array().head( N );
            uStage.array() += dt * A( iStage, i ) * kStages[i].array().tail( N );
        }

        kStages[iStage] = evaluateFunction( rStage, p, uStage );
    }

    r.setZero();
    u.setZero();

    for ( int iStage = 0; iStage < nbStages; iStage++ )
    {
        r.array() += b_RK4[iStage] * kStages[iStage].array().head( N );
        u.array() += b_RK4[iStage] * kStages[iStage].array().tail( N );
    }

    r *= dt;
    u *= dt;
    r += rn;
    u += un;

    // Return area a
    a = r.array() + r0;
    a = boost::math::constants::pi<scalar>() * a.array() * a.array();

    data.col( 0 ) = a;
    this->p = p;
}

fsi::vector TubeFlowExplicitLinearSolidSolver::evaluateFunction(
    const fsi::vector & rn,
    const fsi::vector & pn,
    const fsi::vector & un
    )
{
    assert( rn.rows() == N );
    assert( pn.rows() == N );

    fsi::vector f( 2 * N );
    f.head( N ) = un;

    for ( int i = 0; i < N; i++ )
    {
        if ( i == 0 )
        {
            f( i + N ) = kappa * G / rho * ( rn( i + 1 ) - rn( i ) )
                - E0 / (1 - nu * nu) * rn( i ) / (r0 * r0 * rho)
                + pn( i ) / (rho * h);
        }

        if ( i > 0 && i < N - 1 )
        {
            f( i + N ) = kappa * G / rho * ( rn( i + 1 ) - 2 * rn( i ) + rn( i - 1 ) )
                - E0 / (1 - nu * nu) * rn( i ) / (r0 * r0 * rho)
                + pn( i ) / (rho * h);
        }

        if ( i == N - 1 )
        {
            f( i + N ) = kappa * G / rho * ( -rn( i ) + rn( i - 1 ) )
                - E0 / (1 - nu * nu) * rn( i ) / (r0 * r0 * rho)
                + pn( i ) / (rho * h);
        }
    }

    return f;
}
