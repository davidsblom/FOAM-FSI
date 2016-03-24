
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
        // Ugly hack to get SDC time integration working
        alpha = a0 / (u0 + dx / 0.1);
    }

    SDCTubeFlowFluidSolver::~SDCTubeFlowFluidSolver()
    {}

    scalar SDCTubeFlowFluidSolver::evaluateOutputPressureBoundaryCondition(
        scalar /*pout_n*/,
        scalar /*uout_n*/,
        scalar /*uout*/
        )
    {
        // Ugly hack to get SDC time integration working
        return 0;
    }

    void SDCTubeFlowFluidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & f
        )
    {
        // Check input parameters
        assert( f.rows() == 2 * N );
        assert( a.rows() == N );
        assert( un.rows() == N );
        assert( pn.rows() == N );
        assert( an.rows() == N );
        assert( f.rows() == rhs.rows() );

        f.setZero();

        // Conservation of mass: internal system

        // Left and right face of a: a(i-1/2), a(i+1/2)
        fsi::vector a_lf = 0.5 * ( a.head( N - 2 ) + a.segment( 1, N - 2 ) );
        fsi::vector a_rf = 0.5 * ( a.segment( 1, N - 2 ) + a.segment( 2, N - 2 ) );

        // Left and right face of u: u(i-1/2), u(i+1/2)
        fsi::vector u_lf = 0.5 * ( u.head( N - 2 ) + u.segment( 1, N - 2 ) );
        fsi::vector u_rf = 0.5 * ( u.segment( 1, N - 2 ) + u.segment( 2, N - 2 ) );

        assert( a_lf.rows() == a_rf.rows() );
        assert( u_lf.rows() == u_rf.rows() );
        assert( a_lf.rows() == u_lf.rows() );
        assert( a_lf.rows() == N - 2 );

        // Spatial term of mass conservation
        f.segment( 1, N - 2 ) += u_rf.cwiseProduct( a_rf ) - u_lf.cwiseProduct( a_lf );

        // Add pressure stabilization with coefficient alpha
        f.segment( 1, N - 2 ) -= alpha / rho * ( p.segment( 2, N - 2 ) - 2 * p.segment( 1, N - 2 ) + p.head( N - 2 ) );

        f.segment( 1, N - 2 ) /= dx;

        // Conservation of momentum: internal system

        // Convective term of momentum conservation
        f.segment( N + 1, N - 2 ) += u.segment( 1, N - 2 ).cwiseProduct( u_rf ).cwiseProduct( a_rf ) - u.head( N - 2 ).cwiseProduct( u_lf ).cwiseProduct( a_lf );

        // Pressure term of momentum conservation
        f.segment( N + 1, N - 2 ) += 0.5 * 1.0 / rho * a_rf.cwiseProduct( p.segment( 2, N - 2 ) - p.segment( 1, N - 2 ) );
        f.segment( N + 1, N - 2 ) += 0.5 * 1.0 / rho * a_lf.cwiseProduct( p.segment( 1, N - 2 ) - p.head( N - 2 ) );

        f.segment( N + 1, N - 2 ) /= dx;

        f *= -1.0;
    }

    void SDCTubeFlowFluidSolver::finalizeTimeStep()
    {
        TubeFlowFluidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowFluidSolver::getDOF()
    {
        return 2 * N;
    }

    void SDCTubeFlowFluidSolver::getSolution(
        fsi::vector & solution,
        fsi::vector & /*f*/
        )
    {
        solution.head( N ) = a;
        solution.tail( N ) = a.cwiseProduct( u );

        // Ignore the boundary conditions
        solution( 0 ) = 0;
        solution( N - 1 ) = 0;
        solution( N ) = 0;
        solution( 2 * N - 1 ) = 0;
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

        this->solve( a, p );

        getSolution( result, f );
        evaluateFunction( k, result, t, f );

        finalizeImplicitSolve( k );
    }

    void SDCTubeFlowFluidSolver::prepareImplicitSolve(
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

        pn = pStages.at( kold );
        un = uStages.at( kold );
        an = aStages.at( kold );

        this->rhs = rhs;
    }

    void SDCTubeFlowFluidSolver::finalizeImplicitSolve( int k )
    {
        pStages.at( k + 1 ) = p;
        uStages.at( k + 1 ) = u;
        aStages.at( k + 1 ) = a;
    }

    void SDCTubeFlowFluidSolver::getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        )
    {
        dof.push_back( N );
        enabled.push_back( true );
        names.push_back( "fluid a" );
        dof.push_back( N );
        enabled.push_back( true );
        names.push_back( "fluid a * u" );
    }
}
