
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
        scalar r0
        )
        :
        TubeFlowLinearSolidSolver( N, nu, rho, h, L, dt, G, E0, r0 ),
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
        f.setZero();

        for ( int i = 0; i < N; i++ )
            f( i ) = u( i );

        for ( int i = 1; i < N - 1; i++ )
            f( i + N ) = -1.0 / (rho * h) * ( beta * ( r( i + 1 ) - 2 * r( i ) + r( i - 1 ) ) + gamma * r( i ) - gamma * r0 - p( i ) );
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

        // Ignore the boundary conditions
        solution( N ) = 0;
        solution( 2 * N - 1 ) = 0;
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
        assert( false );
        return 0;
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
            u = uStages.at( k );
            r = rStages.at( k );
        }

        un = uStages.at( kold );
        rn = rStages.at( kold );

        this->rhs = rhs;

        alpha = rho * h / dt;
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
