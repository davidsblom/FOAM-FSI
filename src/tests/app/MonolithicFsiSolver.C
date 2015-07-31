
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MonolithicFsiSolver.H"

namespace tubeflow
{
    MonolithicFsiSolver::MonolithicFsiSolver(
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
        a0( a0 ),
        u0( u0 ),
        p0( p0 ),
        dt( dt ),
        cmk( cmk ),
        N( N ),
        L( L ),
        T( T ),
        rho( rho ),
        un( N ),
        pn( N ),
        an( N ),
        u( N ),
        p( N ),
        a( N ),
        init( false ),
        dx( L / N ),
        alpha( a0 / (u0 + dx / dt) ),
        tau( u0 * dt / L ),
        timeIndex( 0 ),
        t( 0 ),
        p_outn( 0 ),
        p_out( 0 )
    {
        un.fill( u0 );
        an.fill( a0 );
        pn.fill( p0 );
        u.fill( u0 );
        p.fill( p0 );
        a.fill( a0 );
    }

    scalar MonolithicFsiSolver::evaluateInletVelocityBoundaryCondition()
    {
        assert( init );

        return u0 + u0 / 10.0 * std::pow( std::sin( M_PI * timeIndex * tau ), 2 );
    }

    scalar MonolithicFsiSolver::evaluateOutputPressureBoundaryCondition(
        scalar pout_n,
        scalar uout_n,
        scalar uout
        )
    {
        assert( init );

        scalar value = std::sqrt( cmk * cmk - pout_n / (2.0 * rho) );

        value = 2 * rho * ( cmk * cmk - std::pow( value - (uout - uout_n) / 4.0, 2 ) );

        return value;
    }

    void MonolithicFsiSolver::evaluateResidual(
        fsi::vector & x,
        fsi::vector & un,
        fsi::vector & pn,
        fsi::vector & an,
        fsi::vector & R
        )
    {
        // Check input parameters
        assert( init );
        assert( R.rows() == 3 * N );
        assert( x.rows() == 3 * N );
        assert( un.rows() == N );
        assert( pn.rows() == N );
        assert( an.rows() == N );
        assert( t <= T );

        // Determine the velocity and pressure fsi::vectors from the
        // input fsi::vector x
        fsi::vector u = x.head( N );
        fsi::vector p = x.segment( N, N );
        fsi::vector a = x.tail( N );

        // Boundary conditions
        scalar u_in = evaluateInletVelocityBoundaryCondition();
        scalar u_out = 2 * u( N - 2 ) - u( N - 3 );
        scalar u_outn = 2 * un( N - 2 ) - un( N - 3 );
        scalar p_in = 2 * p( 1 ) - p( 2 );
        p_out = evaluateOutputPressureBoundaryCondition( p_outn, u_outn, u_out );

        // Apply boundary conditions

        // Velocity inlet boundary condition inlet
        R( 0 ) = u( 0 ) - u_in;

        // Velocity outlet boundary condition
        R( N - 1 ) = u( N - 1 ) - u_out;

        // Pressure inlet boundary condition
        R( N ) = p( 0 ) - p_in;

        // Pressure outlet boundary condition
        R( 2 * N - 1 ) = p( N - 1 ) - p_out;

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

        // Temporal term of mass conservation
        R.segment( 1, N - 2 ) = dx / dt * ( a.segment( 1, N - 2 ) - an.segment( 1, N - 2 ) );

        // Spatial term of mass conservation
        R.segment( 1, N - 2 ) += u_rf.cwiseProduct( a_rf ) - u_lf.cwiseProduct( a_lf );

        // Add pressure stabilization with coefficient alpha
        R.segment( 1, N - 2 ) -= alpha / rho * ( p.segment( 2, N - 2 ) - 2 * p.segment( 1, N - 2 ) + p.head( N - 2 ) );

        // Conservation of momentum: internal system

        // Temporal term of momentum conservation
        R.segment( N + 1, N - 2 ) = dx / dt * ( u.segment( 1, N - 2 ).cwiseProduct( a.segment( 1, N - 2 ) ) - un.segment( 1, N - 2 ).cwiseProduct( an.segment( 1, N - 2 ) ) );

        // Convective term of momentum conservation
        R.segment( N + 1, N - 2 ) += u.segment( 1, N - 2 ).cwiseProduct( u_rf ).cwiseProduct( a_rf ) - u.head( N - 2 ).cwiseProduct( u_lf ).cwiseProduct( a_lf );

        // Pressure term of momentum conservation
        R.segment( N + 1, N - 2 ) += 0.5 * 1.0 / rho * a_rf.cwiseProduct( p.segment( 2, N - 2 ) - p.segment( 1, N - 2 ) );
        R.segment( N + 1, N - 2 ) += 0.5 * 1.0 / rho * a_lf.cwiseProduct( p.segment( 1, N - 2 ) - p.head( N - 2 ) );

        // Structure part
        for ( int i = 0; i < N; i++ )
        {
            scalar tmp = p0 / (2 * rho) - cmk * cmk;
            tmp /= p( i ) / (2 * rho) - cmk * cmk;
            tmp = a0 * tmp * tmp;

            R( 2 * N + i ) = a( i ) - tmp;
        }
    }

    void MonolithicFsiSolver::finalizeTimeStep()
    {
        assert( init );

        un = u;
        pn = p;
        an = a;
        p_outn = p_out;

        init = false;
    }

    void MonolithicFsiSolver::initTimeStep()
    {
        assert( !init );

        timeIndex++;
        t = timeIndex * dt;

        init = true;
    }

    bool MonolithicFsiSolver::isRunning()
    {
        return t < T;
    }

    void MonolithicFsiSolver::run()
    {
        assert( !init );

        while ( isRunning() )
            solveTimeStep();
    }

    // Generic functor
    template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    struct Functor
    {
        typedef _Scalar Scalar;
        enum {
            InputsAtCompileTime = NX,
            ValuesAtCompileTime = NY
        };
        typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
        typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
        typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

        int m_inputs, m_values;

        Functor()
            :
            m_inputs( InputsAtCompileTime ),
            m_values( ValuesAtCompileTime )
        {}

        virtual ~Functor()
        {}

        Functor(
            int inputs,
            int values
            ) : m_inputs( inputs ), m_values( values ) {}

        int inputs() const
        {
            return m_inputs;
        }

        int values() const
        {
            return m_values;
        }
    };

    struct residualFunctor : Functor<scalar>
    {
        residualFunctor(
            MonolithicFsiSolver * fsi,
            fsi::vector * un,
            fsi::vector * pn,
            fsi::vector * an
            )
            :
            Functor<scalar>( 3 * fsi->N, 3 * fsi->N ),
            fsi( fsi ),
            un( un ),
            pn( pn ),
            an( an )
        {}

        virtual ~residualFunctor()
        {}

        int operator()(
            fsi::vector & x,
            fsi::vector & fvec
            ) const
        {
            fsi->evaluateResidual( x, *un, *pn, *an, fvec );
            return 0;
        }

        MonolithicFsiSolver * fsi;
        fsi::vector * un;
        fsi::vector * pn;
        fsi::vector * an;
    };

    void MonolithicFsiSolver::solve()
    {
        assert( init );

        // Find the velocity, pressure and displacement for which the residual is zero

        // Initial guess is current solution or the solution at previous time step
        fsi::vector x( 3 * N );
        x.head( N ) = u;
        x.segment( N, N ) = p;
        x.tail( N ) = a;

        // Optimize the residual function
        residualFunctor functor( this, &un, &pn, &an );
        Eigen::NumericalDiff<residualFunctor, Eigen::Central> numDiff( functor );
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<residualFunctor, Eigen::Central>, scalar> lm( numDiff );

        lm.parameters.maxfev = 2000;
        lm.parameters.xtol = 1.0e-13;
        lm.parameters.ftol = 1.0e-13;

        int ret = lm.minimize( x );

        // 2: RelativeErrorTooSmall
        assert( ret == 2 );

        // Save the velocity, pressure and displacement
        u = x.head( N );
        p = x.segment( N, N );
        a = x.tail( N );
    }

    void MonolithicFsiSolver::solveTimeStep()
    {
        assert( !init );

        initTimeStep();
        solve();
        finalizeTimeStep();
    }
}
