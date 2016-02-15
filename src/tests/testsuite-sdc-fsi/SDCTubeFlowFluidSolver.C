
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCTubeFlowFluidSolver.H"

namespace tubeflow
{
    SDCTubeFlowFluidSolver::SDCTubeFlowFluidSolver(
        scalar a0,
        scalar u0,
        scalar p0,
        scalar dt,
        scalar cmk,
        int N,
        scalar L,
        scalar T,
        scalar rho
    )
    :
    TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho, true ),
    k( 0 ),
    pStages(),
    uStages(),
    aStages()
    {

    }

    SDCTubeFlowFluidSolver::~SDCTubeFlowFluidSolver()
    {

    }

    void SDCTubeFlowFluidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & f
        )
    {
        f.setZero();
    }

    void SDCTubeFlowFluidSolver::finalizeTimeStep()
    {
        TubeFlowFluidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowFluidSolver::getDOF()
    {
        return 2*N;
    }

    void SDCTubeFlowFluidSolver::getSolution(
        fsi::vector & solution,
        fsi::vector & /*f*/
        )
    {
        solution.head( N ) = a;
        solution.tail( N ) = a.cwiseProduct( u );
    }

    void SDCTubeFlowFluidSolver::setSolution(
        const fsi::vector & /*solution*/,
        const fsi::vector & /*f*/
        )
    {
        assert( false );
    }

    scalar SDCTubeFlowFluidSolver::getEndTime()
    {
        return T;
    }

    scalar SDCTubeFlowFluidSolver::getTimeStep()
    {
        return dt;
    }

    void SDCTubeFlowFluidSolver::nextTimeStep()
    {
        for ( int i = 0; i < k; i++ )
        {
            pStages.at( i ) = p;
            uStages.at( i ) = u;
            aStages.at( i ) = a;
        }
    }

    void SDCTubeFlowFluidSolver::initTimeStep()
    {
        TubeFlowFluidSolver::initTimeStep();
    }

    void SDCTubeFlowFluidSolver::setNumberOfImplicitStages( int k )
    {
        this->k = k + 1;

        pStages.clear();
        uStages.clear();
        aStages.clear();

        for ( int i = 0; i < k + 1; i++ )
        {
            pStages.push_back( p );
            uStages.push_back( u );
            aStages.push_back( a );
        }
    }

    void SDCTubeFlowFluidSolver::implicitSolve(
        bool corrector,
        const int k,
        const int /*kold*/,
        const scalar t,
        const scalar dt,
        const fsi::vector & /*qold*/,
        const fsi::vector & rhs,
        fsi::vector & f,
        fsi::vector & result
    )
    {
        this->dt = dt;
        this->t = t;
        this->tau = u0 * dt / L;

        if ( corrector )
        {
            p = pStages.at( k + 1 );
            u = uStages.at( k + 1 );
            a = aStages.at( k + 1 );
        }

        if ( not corrector )
        {
            p = pStages.at( k );
            u = uStages.at( k );
            a = aStages.at( k );
        }

        pn = pStages.at( k );
        un = uStages.at( k );
        an = aStages.at( k );

        this->rhs = rhs;

        this->solve( a, p );

        pStages.at( k + 1 ) = p;
        uStages.at( k + 1 ) = u;
        aStages.at( k + 1 ) = a;

        getSolution( result, f );
    }
}
