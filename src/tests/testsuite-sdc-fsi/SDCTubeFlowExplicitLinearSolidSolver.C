
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCTubeFlowExplicitLinearSolidSolver.H"

namespace tubeflow
{
    SDCTubeFlowExplicitLinearSolidSolver::SDCTubeFlowExplicitLinearSolidSolver(
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
        scalar T,
        std::shared_ptr<sdc::DataStorage> dataStorage
        )
        :
        TubeFlowExplicitLinearSolidSolver( N, nu, rho, h, L, dt, G, E0, r0, p0, T, dataStorage ),
        k( 0 )
    {}

    SDCTubeFlowExplicitLinearSolidSolver::SDCTubeFlowExplicitLinearSolidSolver(
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
        TubeFlowExplicitLinearSolidSolver( N, nu, rho, h, L, dt, G, E0, r0, p0, T ),
        k( 0 )
    {}

    SDCTubeFlowExplicitLinearSolidSolver::~SDCTubeFlowExplicitLinearSolidSolver()
    {}

    void SDCTubeFlowExplicitLinearSolidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & /*f*/
        )
    {}

    void SDCTubeFlowExplicitLinearSolidSolver::finalizeTimeStep()
    {
        TubeFlowExplicitLinearSolidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowExplicitLinearSolidSolver::getDOF()
    {
        return 0;
    }

    void SDCTubeFlowExplicitLinearSolidSolver::getSolution(
        fsi::vector & /*solution*/,
        fsi::vector & /*f*/
        )
    {}

    void SDCTubeFlowExplicitLinearSolidSolver::setSolution(
        const fsi::vector & /*solution*/,
        const fsi::vector &     /*f*/
        )
    {
        assert( false );
    }

    scalar SDCTubeFlowExplicitLinearSolidSolver::getEndTime()
    {
        return T;
    }

    scalar SDCTubeFlowExplicitLinearSolidSolver::getTimeStep()
    {
        return dt;
    }

    void SDCTubeFlowExplicitLinearSolidSolver::nextTimeStep()
    {
        for ( int i = 0; i < k; i++ )
        {
            uStages.at( i ) = u;
            rStages.at( i ) = r;
            pStages.at( i ) = p;
        }
    }

    void SDCTubeFlowExplicitLinearSolidSolver::initTimeStep()
    {
        TubeFlowExplicitLinearSolidSolver::initTimeStep();
    }

    void SDCTubeFlowExplicitLinearSolidSolver::setNumberOfImplicitStages( int k )
    {
        this->k = k + 1;

        uStages.clear();
        rStages.clear();
        pStages.clear();

        for ( int i = 0; i < k + 1; i++ )
        {
            uStages.push_back( u );
            rStages.push_back( r );
            pStages.push_back( p );
        }
    }

    void SDCTubeFlowExplicitLinearSolidSolver::implicitSolve(
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

    void SDCTubeFlowExplicitLinearSolidSolver::prepareImplicitSolve(
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
        TubeFlowExplicitLinearSolidSolver::k = k;
        TubeFlowExplicitLinearSolidSolver::corrector = corrector;

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
    }

    void SDCTubeFlowExplicitLinearSolidSolver::finalizeImplicitSolve( int k )
    {
        uStages.at( k + 1 ) = u;
        rStages.at( k + 1 ) = r;
        pStages.at( k + 1 ) = p;
    }

    void SDCTubeFlowExplicitLinearSolidSolver::getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        )
    {
        dof.push_back( 0 );
        enabled.push_back( false );
        names.push_back( "solid" );
    }
}
