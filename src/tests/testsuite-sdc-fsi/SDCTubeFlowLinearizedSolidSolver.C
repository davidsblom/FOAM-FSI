
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCTubeFlowLinearizedSolidSolver.H"

namespace tubeflow
{
    SDCTubeFlowLinearizedSolidSolver::SDCTubeFlowLinearizedSolidSolver(
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
        TubeFlowLinearizedSolidSolver( N, nu, rho, h, L, dt, G, E0, r0 )
    {}

    SDCTubeFlowLinearizedSolidSolver::~SDCTubeFlowLinearizedSolidSolver()
    {}

    void SDCTubeFlowLinearizedSolidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & f
        )
    {
        f.setZero();

        for ( int i = 0; i < N; i++ )
            f( i ) = u( i );

        for ( int i = 1; i < N - 1; i++ )
            f( i + N ) = -1.0 / (rho * h) * ( beta * ( r( i + 1 ) - 2 * r( i ) + r( i - 1 ) ) + gamma * r( i ) - p( i ) );
    }

    void SDCTubeFlowLinearizedSolidSolver::finalizeTimeStep()
    {
        TubeFlowLinearizedSolidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowLinearizedSolidSolver::getDOF()
    {
        return 2 * N;
    }

    void SDCTubeFlowLinearizedSolidSolver::getSolution(
        fsi::vector & solution,
        fsi::vector &     /*f*/
        )
    {
        solution.head( N ) = r;
        solution.tail( N ) = u;

        // Ignore the boundary conditions
        solution( N ) = 0;
        solution( 2 * N - 1 ) = 0;
    }

    void SDCTubeFlowLinearizedSolidSolver::setSolution(
        const fsi::vector & /*solution*/,
        const fsi::vector &     /*f*/
        )
    {
        assert( false );
    }

    scalar SDCTubeFlowLinearizedSolidSolver::getEndTime()
    {
        assert( false );
        return 0;
    }

    scalar SDCTubeFlowLinearizedSolidSolver::getTimeStep()
    {
        return dt;
    }

    void SDCTubeFlowLinearizedSolidSolver::nextTimeStep()
    {
        for ( int i = 0; i < k; i++ )
        {
            uStages.at( i ) = u;
            rStages.at( i ) = r;
        }
    }

    void SDCTubeFlowLinearizedSolidSolver::initTimeStep()
    {
        TubeFlowLinearizedSolidSolver::initTimeStep();
    }

    void SDCTubeFlowLinearizedSolidSolver::setNumberOfImplicitStages( int k )
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

    void SDCTubeFlowLinearizedSolidSolver::implicitSolve(
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

        getSolution( result, f );
        evaluateFunction( k, result, t, f );

        finalizeImplicitSolve( k );
    }

    void SDCTubeFlowLinearizedSolidSolver::prepareImplicitSolve(
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
            u = uStages.at( k );
            r = rStages.at( k );
        }

        un = uStages.at( kold );
        rn = rStages.at( kold );

        this->rhs = rhs;

        alpha = rho * h / dt;
        factorizeMatrix();
    }

    void SDCTubeFlowLinearizedSolidSolver::finalizeImplicitSolve( int k )
    {
        uStages.at( k + 1 ) = u;
        rStages.at( k + 1 ) = r;
    }

    void SDCTubeFlowLinearizedSolidSolver::getVariablesInfo(
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
