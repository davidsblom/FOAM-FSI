
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCTubeFlowLinearSolidSolver.H"

namespace tubeflow
{
    SDCTubeFlowLinearSolidSolver::SDCTubeFlowLinearSolidSolver(
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
        TubeFlowLinearSolidSolver( N, nu, rho, h, L, dt, G, E0, r0, T ),
        k( 0 ),
        uStages(),
        rStages()
    {}

    SDCTubeFlowLinearSolidSolver::~SDCTubeFlowLinearSolidSolver()
    {}

    void SDCTubeFlowLinearSolidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & f
        )
    {
        for ( int i = 0; i < N; i++ )
            f( i ) = u( i );

        for ( int i = 1; i < N - 1; i++ )
            f( i + N ) = kappa * G / (dx * dx * rho) * ( r( i + 1 ) - 2 * r( i ) + r( i - 1 ) )
                - E0 / (1 - nu * nu) * r( i ) / (r0 * r0 * rho)
                + p( i ) / (rho * h)
                + E0 / (rho * r0) / (1 - nu * nu);

        f( N ) = kappa * G / (dx * dx * rho) * ( r( 1 ) - r( 0 ) )
            - E0 / (1 - nu * nu) * r( 0 ) / (r0 * r0 * rho)
            + p( 0 ) / (rho * h)
            + E0 / (rho * r0) / (1 - nu * nu);
        f( 2 * N - 1 ) = kappa * G / (dx * dx * rho) * ( -r( N - 1 ) + r( N - 2 ) )
            - E0 / (1 - nu * nu) * r( N - 1 ) / (r0 * r0 * rho)
            + p( N - 1 ) / (rho * h)
            + E0 / (rho * r0) / (1 - nu * nu);
    }

    void SDCTubeFlowLinearSolidSolver::finalizeTimeStep()
    {
        TubeFlowLinearSolidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowLinearSolidSolver::getDOF()
    {
        return 2 * N;
    }

    void SDCTubeFlowLinearSolidSolver::getSolution(
        fsi::vector & solution,
        fsi::vector & /*f*/
        )
    {
        solution.head( N ) = r;
        solution.tail( N ) = u;
    }

    void SDCTubeFlowLinearSolidSolver::setSolution(
        const fsi::vector & /*solution*/,
        const fsi::vector & /*f*/
        )
    {
        assert( false );
    }

    scalar SDCTubeFlowLinearSolidSolver::getEndTime()
    {
        return T;
    }

    scalar SDCTubeFlowLinearSolidSolver::getTimeStep()
    {
        return dt;
    }

    void SDCTubeFlowLinearSolidSolver::nextTimeStep()
    {
        for ( int i = 0; i < k; i++ )
        {
            uStages.at( i ) = u;
            rStages.at( i ) = r;
        }
    }

    void SDCTubeFlowLinearSolidSolver::initTimeStep()
    {
        TubeFlowLinearSolidSolver::initTimeStep();
    }

    void SDCTubeFlowLinearSolidSolver::setNumberOfImplicitStages( int k )
    {
        this->k = k + 1;

        uStages.clear();
        rStages.clear();

        for ( int i = 0; i < k + 1; i++ )
        {
            uStages.push_back( u );
            rStages.push_back( r );
        }
    }

    void SDCTubeFlowLinearSolidSolver::implicitSolve(
        bool corrector,
        const int k,
        const int kold,
        const scalar t,
        const scalar dt,
        const fsi::vector & qold,
        const fsi::vector & rhs,
        fsi::vector & f,
        fsi::vector & result
        )
    {
        prepareImplicitSolve( corrector, k, kold, t, dt, qold, rhs );

        fsi::vector a;
        this->solve( p, a );

        getSolution( result, f );
        evaluateFunction( k, result, t, f );

        finalizeImplicitSolve( k );
    }

    void SDCTubeFlowLinearSolidSolver::prepareImplicitSolve(
        bool corrector,
        const int k,
        const int kold,
        const scalar t,
        const scalar dt,
        const fsi::vector & /*qold*/,
        const fsi::vector & rhs
        )
    {
        this->dt = dt;
        this->t = t;

        if ( corrector )
        {
            u = uStages.at( k + 1 );
            r = rStages.at( k + 1 );
        }

        if ( not corrector )
        {
            u = uStages.at( kold );
            r = rStages.at( kold );
        }

        un = uStages.at( kold );
        rn = rStages.at( kold );

        this->rhs = rhs;

        factorizeMatrix();
    }

    void SDCTubeFlowLinearSolidSolver::finalizeImplicitSolve( int k )
    {
        uStages.at( k + 1 ) = u;
        rStages.at( k + 1 ) = r;
    }

    void SDCTubeFlowLinearSolidSolver::getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        )
    {
        dof.push_back( N );
        enabled.push_back( true );
        names.push_back( "solid u" );
        dof.push_back( N );
        enabled.push_back( true );
        names.push_back( "solid r" );
    }
}
