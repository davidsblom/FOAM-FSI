
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCTubeFlowBDFLinearizedSolidSolver.H"

namespace tubeflow
{
    SDCTubeFlowBDFLinearizedSolidSolver::SDCTubeFlowBDFLinearizedSolidSolver(
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
        TubeFlowLinearizedSolidSolver( N, nu, rho, h, L, dt, G, E0, r0, T, timeOrder )
    {}

    SDCTubeFlowBDFLinearizedSolidSolver::~SDCTubeFlowBDFLinearizedSolidSolver()
    {}

    void SDCTubeFlowBDFLinearizedSolidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & /*f*/
        )
    {}

    void SDCTubeFlowBDFLinearizedSolidSolver::finalizeTimeStep()
    {
        unStages = uStages;
        rnStages = rStages;
        TubeFlowLinearizedSolidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowBDFLinearizedSolidSolver::getDOF()
    {
        return 0;
    }

    void SDCTubeFlowBDFLinearizedSolidSolver::getSolution(
        fsi::vector & /*solution*/,
        fsi::vector &     /*f*/
        )
    {}

    void SDCTubeFlowBDFLinearizedSolidSolver::setSolution(
        const fsi::vector & /*solution*/,
        const fsi::vector &     /*f*/
        )
    {
        assert( false );
    }

    scalar SDCTubeFlowBDFLinearizedSolidSolver::getEndTime()
    {
        return T;
    }

    scalar SDCTubeFlowBDFLinearizedSolidSolver::getTimeStep()
    {
        return dt;
    }

    void SDCTubeFlowBDFLinearizedSolidSolver::nextTimeStep()
    {
        for ( int i = 0; i < k; i++ )
        {
            uStages.at( i ) = u;
            rStages.at( i ) = r;
        }
    }

    void SDCTubeFlowBDFLinearizedSolidSolver::initTimeStep()
    {
        TubeFlowLinearizedSolidSolver::initTimeStep();
    }

    void SDCTubeFlowBDFLinearizedSolidSolver::setNumberOfImplicitStages( int k )
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

    void SDCTubeFlowBDFLinearizedSolidSolver::implicitSolve(
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

    void SDCTubeFlowBDFLinearizedSolidSolver::prepareImplicitSolve(
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

        TubeFlowLinearizedSolidSolver::rStages[0] = rn;
        TubeFlowLinearizedSolidSolver::uStages[0] = un;

        if ( timeIndex == 1 && timeOrder == 2 && k > 1 )
        {
            TubeFlowLinearizedSolidSolver::rStages.clear();
            TubeFlowLinearizedSolidSolver::uStages.clear();

            TubeFlowLinearizedSolidSolver::rStages.push_back( rStages.at( kold - 1 ) );
            TubeFlowLinearizedSolidSolver::uStages.push_back( uStages.at( kold - 1 ) );
            TubeFlowLinearizedSolidSolver::rStages.push_back( rn );
            TubeFlowLinearizedSolidSolver::uStages.push_back( un );
            TubeFlowLinearizedSolidSolver::rStages.push_back( r );
            TubeFlowLinearizedSolidSolver::uStages.push_back( u );

            alpha = {
                1. / 3., -4. / 3., 1.
            };
            beta = 2. / 3.;
        }

        if ( timeIndex > 1 && timeOrder == 2 )
        {
            if ( kold - 1 < 0 )
            {
                TubeFlowLinearizedSolidSolver::rStages[0] = rnStages.at( rnStages.size() - 2 );
                TubeFlowLinearizedSolidSolver::uStages[0] = unStages.at( unStages.size() - 2 );
            }
            else
            {
                TubeFlowLinearizedSolidSolver::rStages[0] = rStages.at( kold - 1 );
                TubeFlowLinearizedSolidSolver::uStages[0] = uStages.at( kold - 1 );
            }

            TubeFlowLinearizedSolidSolver::rStages[1] = rn;
            TubeFlowLinearizedSolidSolver::uStages[1] = un;
        }

        this->rhs = rhs;

        factorizeMatrix();
    }

    void SDCTubeFlowBDFLinearizedSolidSolver::finalizeImplicitSolve( int k )
    {
        uStages.at( k + 1 ) = u;
        rStages.at( k + 1 ) = r;
    }

    void SDCTubeFlowBDFLinearizedSolidSolver::getVariablesInfo(
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
