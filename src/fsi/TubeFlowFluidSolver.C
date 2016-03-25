
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowFluidSolver.H"
#include <iomanip>

namespace tubeflow
{
    TubeFlowFluidSolver::TubeFlowFluidSolver(
        scalar a0,
        scalar u0,
        scalar p0,
        scalar dt,
        scalar cmk,
        int N,
        scalar L,
        scalar T,
        scalar rho
        ) :
        BaseMultiLevelSolver( N, 1, p0 ),
        a0( a0 ),
        u0( u0 ),
        p0( p0 ),
        dt( dt ),
        dx( L / N ),
        cmk( cmk ),
        rho( rho ),
        L( L ),
        T( T ),
        alpha( a0 / (u0 + dx / dt) ),
        tau( u0 * dt / L ),
        p_outn( 0 ),
        p_out( 0 ),
        un( N ),
        pn( N ),
        an( N ),
        u( N ),
        p( N ),
        a( N ),
        rhs( 2 * N ),
        iter( 0 ),
        minIter( 1 ),
        maxIter( 30 ),
        tol( 1.0e-14 ),
        nbRes( 0 ),
        nbJac( 0 ),
        grid(),
        diffJacobian( false )
    {
        // Verify input parameters
        assert( dt > 0 );
        assert( a0 > 0 );
        assert( cmk > 0 );
        assert( N > 0 );
        assert( L > 0 );
        assert( T > 0 );
        assert( rho > 0 );
        assert( alpha > 0 );
        assert( tau > 0 );

        un.fill( u0 );
        an.fill( a0 );
        pn.fill( p0 );
        u.fill( u0 );
        p.fill( p0 );
        a.fill( a0 );
        data.fill( p0 );
        rhs.setZero();
    }

    TubeFlowFluidSolver::TubeFlowFluidSolver(
        scalar a0,
        scalar u0,
        scalar p0,
        scalar dt,
        scalar cmk,
        int N,
        scalar L,
        scalar T,
        scalar rho,
        bool diffJacobian
        )
        :
        BaseMultiLevelSolver( N, 1, p0 ),
        a0( a0 ),
        u0( u0 ),
        p0( p0 ),
        dt( dt ),
        dx( L / N ),
        cmk( cmk ),
        rho( rho ),
        L( L ),
        T( T ),
        alpha( a0 / (u0 + dx / dt) ),
        tau( u0 * dt / L ),
        p_outn( 0 ),
        p_out( 0 ),
        un( N ),
        pn( N ),
        an( N ),
        u( N ),
        p( N ),
        a( N ),
        rhs( 2 * N ),
        iter( 0 ),
        minIter( 1 ),
        maxIter( 30 ),
        tol( 1.0e-14 ),
        nbRes( 0 ),
        nbJac( 0 ),
        grid(),
        diffJacobian( diffJacobian )
    {
        // Verify input parameters
        assert( dt > 0 );
        assert( a0 > 0 );
        assert( cmk > 0 );
        assert( N > 0 );
        assert( L > 0 );
        assert( T > 0 );
        assert( rho > 0 );
        assert( alpha > 0 );
        assert( tau > 0 );

        un.fill( u0 );
        an.fill( a0 );
        pn.fill( p0 );
        u.fill( u0 );
        p.fill( p0 );
        a.fill( a0 );
        data.fill( p0 );
        rhs.setZero();
    }

    TubeFlowFluidSolver::~TubeFlowFluidSolver()
    {}

    void TubeFlowFluidSolver::calcGrid()
    {
        if ( grid.rows() != N )
        {
            grid.resize( N, 1 );
            grid.setZero();

            for ( int i = 0; i < N; i++ )
                grid( i, 0 ) = dx * i + 0.5 * dx;
        }
    }

    scalar TubeFlowFluidSolver::evaluateInletVelocityBoundaryCondition()
    {
        return u0 + u0 / 10.0 * std::pow( std::sin( M_PI * t * u0 / L ), 2 );
    }

    scalar TubeFlowFluidSolver::evaluateOutputPressureBoundaryCondition(
        scalar pout_n,
        scalar uout_n,
        scalar uout
        )
    {
        scalar value = std::sqrt( cmk * cmk - pout_n / (2.0 * rho) );

        value = 2 * rho * ( cmk * cmk - std::pow( value - (uout - uout_n) / 4.0, 2 ) );

        return value;
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

        Functor(
            int inputs,
            int values
            )
            :
            m_inputs( inputs ),
            m_values( values )
        {}

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
            TubeFlowFluidSolver * fluid,
            const fsi::vector * a,
            const fsi::vector * un,
            const fsi::vector * pn,
            const fsi::vector * an
            )
            :
            Functor<scalar>( 2 * fluid->N, 2 * fluid->N )
        {
            this->fluid = fluid;
            this->a = a;
            this->un = un;
            this->pn = pn;
            this->an = an;
        }

        int operator()(
            fsi::vector & x,
            fsi::vector & fvec
            ) const
        {
            fluid->evaluateResidual( x, *a, *un, *pn, *an, fvec );
            return 0;
        }

        TubeFlowFluidSolver * fluid;
        const fsi::vector * a;
        const fsi::vector * un;
        const fsi::vector * pn;
        const fsi::vector * an;
    };

    void TubeFlowFluidSolver::evaluateJacobian(
        const fsi::vector & x,
        const fsi::vector & a,
        const fsi::vector & un,
        const fsi::vector & pn,
        const fsi::vector & an,
        matrix & J
        )
    {
        if ( diffJacobian )
        {
            residualFunctor functor( this, &a, &un, &pn, &an );
            Eigen::NumericalDiff<residualFunctor, Eigen::Central> numDiff( functor );

            numDiff.df( x, J );
            return;
        }

        // Check input parameters
        assert( u.rows() == a.rows() );
        assert( x.rows() == u.rows() * 2 );
        assert( x.rows() == a.rows() * 2 );
        assert( p.rows() == u.rows() );
        assert( u.rows() == N );
        assert( x.rows() == 2 * N );
        assert( J.rows() == J.cols() );
        assert( J.rows() == 2 * N );
        assert( init );
        nbJac++;

        // Set Jacobian to zero since the Jacobian is a sparse matrix
        J.setZero();

        // Determine the velocity and pressure fsi::vectors from
        // the input fsi::vector x
        fsi::vector u = x.head( N );
        fsi::vector p = x.tail( N );

        for ( int i = 0; i < N; i++ )
        {
            if ( i == 0 )
            {
                // Derivatives mass conservation
                J( 0, 0 ) = 1.0;

                J( N, N ) = 1.0;
                J( N, N + 1 ) = -2.0;
                J( N, N + 2 ) = 1.0;
            }

            if ( i == N - 1 )
            {
                // Derivatives mass conservation
                J( N - 1, N - 3 ) = 1.0;
                J( N - 1, N - 2 ) = -2.0;
                J( N - 1, N - 1 ) = 1.0;

                J( 2 * N - 1, N - 3 ) = 0.5 * std::sqrt( 2 ) * std::sqrt( 2 * cmk * cmk * rho - p_outn ) * std::sqrt( rho );
                J( 2 * N - 1, N - 3 ) += -0.5 * u( N - 2 ) * rho;
                J( 2 * N - 1, N - 3 ) += 0.25 * u( N - 3 ) * rho;
                J( 2 * N - 1, N - 3 ) += 0.5 * un( N - 2 ) * rho;
                J( 2 * N - 1, N - 3 ) += -0.25 * un( N - 3 ) * rho;

                J( 2 * N - 1, N - 2 ) = -std::sqrt( 2 ) * std::sqrt( 2 * cmk * cmk * rho - p_outn ) * std::sqrt( rho );
                J( 2 * N - 1, N - 2 ) += u( N - 2 ) * rho;
                J( 2 * N - 1, N - 2 ) += -0.5 * u( N - 3 ) * rho;
                J( 2 * N - 1, N - 2 ) += -un( N - 2 ) * rho;
                J( 2 * N - 1, N - 2 ) += 0.5 * un( N - 3 ) * rho;

                J( 2 * N - 1, 2 * N - 1 ) = 1.0;
            }

            if ( i > 0 && i < N - 1 )
            {
                // Derivatives mass conservation
                J( i, i - 1 ) = 0.25 * ( -a( i - 1 ) - a( i ) );
                J( i, i ) = 0.25 * ( -a( i - 1 ) + a( i + 1 ) );
                J( i, i + 1 ) = 0.25 * ( a( i ) + a( i + 1 ) );
                J( i, N + i - 1 ) = -alpha / rho;
                J( i, N + i ) = 2.0 * alpha / rho;
                J( i, N + i + 1 ) = -alpha / rho;

                // Derivatives momentum conservation

                J( N + i, i - 1 ) = -0.5 * u( i - 1 ) * a( i - 1 );
                J( N + i, i - 1 ) += -0.5 * u( i - 1 ) * a( i );
                J( N + i, i - 1 ) += -0.25 * u( i ) * a( i - 1 );
                J( N + i, i - 1 ) += -0.25 * u( i ) * a( i );

                J( N + i, i ) = dx / dt * a( i );
                J( N + i, i ) += 0.5 * u( i ) * a( i );
                J( N + i, i ) += 0.25 * u( i + 1 ) * a( i );
                J( N + i, i ) += 0.5 * u( i ) * a( i + 1 );
                J( N + i, i ) += 0.25 * u( i + 1 ) * a( i + 1 );
                J( N + i, i ) += -0.25 * u( i - 1 ) * a( i );
                J( N + i, i ) += -0.25 * u( i - 1 ) * a( i - 1 );

                J( N + i, i + 1 ) = 0.25 * u( i ) * a( i ) + 0.25 * u( i ) * a( i + 1 );

                J( N + i, N + i - 1 ) = 1.0 / rho * ( -0.25 * a( i - 1 ) - 0.25 * a( i ) );
                J( N + i, N + i ) = 1.0 / rho * ( -0.25 * a( i + 1 ) + 0.25 * a( i - 1 ) );
                J( N + i, N + i + 1 ) = 1.0 / rho * ( 0.25 * a( i ) + 0.25 * a( i + 1 ) );
            }
        }
    }

    /*
     * Evaluate the residual function.
     * Input parameters: solution fsi::vector x, solution fsi::vector a
     * Output: residual fsi::vector R
     */
    void TubeFlowFluidSolver::evaluateResidual(
        const fsi::vector & x,
        const fsi::vector & a,
        const fsi::vector & un,
        const fsi::vector & /*pn*/,
        const fsi::vector & an,
        fsi::vector & R
        )
    {
        // Check input parameters
        assert( init );
        assert( R.rows() == 2 * N );
        assert( x.rows() == 2 * N );
        assert( a.rows() == N );
        assert( un.rows() == N );
        assert( an.rows() == N );
        assert( R.rows() == rhs.rows() );

        R.setZero();
        nbRes++;

        // Determine the velocity and pressure fsi::vectors from the
        // input fsi::vector x
        fsi::vector u = x.head( N );
        fsi::vector p = x.tail( N );

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

        // SDC terms
        R.segment( 1, N - 2 ) -= (dx / dt) * rhs.segment( 1, N - 2 );
        R.segment( N + 1, N - 2 ) -= (dx / dt) * rhs.segment( N + 1, N - 2 );
    }

    void TubeFlowFluidSolver::finalizeTimeStep()
    {
        assert( init );

        un = u;
        pn = p;
        an = a;
        p_outn = p_out;

        init = false;
    }

    void TubeFlowFluidSolver::getWritePositions( matrix & writePositions )
    {
        calcGrid();

        writePositions = grid;
    }

    void TubeFlowFluidSolver::getReadPositions( matrix & readPositions )
    {
        calcGrid();

        readPositions = grid;
    }

    void TubeFlowFluidSolver::initTimeStep()
    {
        assert( !init );

        timeIndex++;
        t = timeIndex * dt;

        init = true;
    }

    bool TubeFlowFluidSolver::isConvergence( const fsi::vector & R )
    {
        assert( init );
        assert( iter <= maxIter );

        scalar norm = R.norm();

        if ( std::isinf( norm ) )
        {
            std::string msg = "The residual function of the fluid solver contains infinite numbers. Unable to continue the computation.";
            std::cout << msg << std::endl;
            throw std::string( msg );
        }

        if ( std::isnan( norm ) )
        {
            std::string msg = "The residual function of the fluid solver contains NaNs. Unable to continue the computation.";
            std::cout << msg << std::endl;
            throw std::string( msg );
        }

        if ( iter < minIter )
            return false;

        if ( iter >= maxIter )
            return true;

        if ( R.cwiseAbs().maxCoeff() > tol )
            return false;

        return norm < tol;
    }

    bool TubeFlowFluidSolver::isRunning()
    {
        assert( !init );

        return t < T;
    }

    void TubeFlowFluidSolver::resetSolution()
    {}

    void TubeFlowFluidSolver::solve(
        const matrix & input,
        matrix & output
        )
    {
        assert( input.rows() == N );
        assert( input.cols() == 1 );

        // Map the matrices to fsi::vectors
        fsi::vector a = Eigen::Map<const fsi::vector> ( input.data(), input.rows() );
        fsi::vector p;

        // Call solve with the correct arguments
        solve( a, p );

        output = p;

        assert( output.rows() == N );
        assert( output.cols() == 1 );
    }

    void TubeFlowFluidSolver::solve(
        const fsi::vector & a,
        fsi::vector & p
        )
    {
        assert( init );
        assert( a.rows() == N );

        std::cout << "Solving fluid domain with size " << N << std::endl;

        // Find the velocity and pressure for which the residual is zero

        // Initial guess is current solution or the solution at previous time step
        fsi::vector x( 2 * N );
        x.head( N ) = u;
        x.tail( N ) = this->p;

        // Optimize the residual function

        iter = 0;

        // Wolfe conditions settings
        scalar c1 = 1.0e-4;
        scalar c2 = 0.9;
        int nbBackTrack = 0;

        // Initialize variables
        fsi::vector R( 2 * N ), Rtry( 2 * N ), dx( 2 * N );
        matrix J( 2 * N, 2 * N ), Jtry( 2 * N, 2 * N );

        // Evaluate the residual
        evaluateResidual( x, a, un, pn, an, R );

        // Newton-Raphson solver
        while ( !isConvergence( R ) )
        {
            iter++;

            evaluateJacobian( x, a, un, pn, an, J );

            dx = J.partialPivLu().solve( -R );

            evaluateResidual( x + dx, a, un, pn, an, Rtry );

            if ( isConvergence( Rtry ) )
            {
                R = Rtry;
                x += dx;
                break;
            }

            // Evaluate Armijo rule and strong Wolfe condition

            bool armijo = false;
            bool wolfe = false;
            scalar alpha = 1.0;

            while ( !armijo || !wolfe )
            {
                evaluateResidual( x.array() + alpha * dx.array(), a, un, pn, an, Rtry );

                armijo = Rtry.norm() <= ( R + (c1 * alpha * dx.transpose() * J).transpose() ).norm();

                evaluateJacobian( x.array() + alpha * dx.array(), a, un, pn, an, Jtry );

                wolfe = (Jtry * dx).norm() >= c2 * (J * dx).norm();

                if ( !armijo || !wolfe )
                {
                    alpha *= 0.9;
                    nbBackTrack++;
                }

                // Fall back to Levenberg-Marquardt algorithm when exact Newton method is not converging
                if ( alpha < 0.0001 )
                {
                    std::string msg = "Newton-Raphson method is not converging. Unable to solve fluid problem.";
                    std::cout << msg << std::endl;

                    // Use the Levenberg Marquardt algorithm to solve the non-linear problem.
                    residualFunctor functor( this, &a, &un, &pn, &an );
                    Eigen::NumericalDiff<residualFunctor, Eigen::Central> numDiff( functor );
                    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<residualFunctor, Eigen::Central>, scalar> lm( numDiff );

                    lm.parameters.maxfev = 2000;
                    lm.parameters.xtol = 1.0e-13;
                    lm.parameters.ftol = 1.0e-13;

                    fsi::vector x0 = x;
                    int ret = lm.minimize( x );
                    alpha = 1;
                    dx = x0 - x;

                    // 2: RelativeErrorTooSmall
                    assert( ret == 2 );

                    if ( ret != 2 )
                        throw std::runtime_error( "The Levenberg Marquardt solver has not converged." );

                    // throw std::runtime_error( msg );
                }
            }

            x += alpha * dx;

            evaluateResidual( x, a, un, pn, an, R );
        }

        assert( isConvergence( R ) );

        // Save the velocity and the pressure
        u = x.head( N );
        this->p = x.tail( N );
        this->a = a;
        p = this->p;
        data.col( 0 ) = p;
    }
}
