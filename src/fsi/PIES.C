
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PIES.H"
#include <math.h>

namespace sdc
{
    PIES::PIES(
        std::shared_ptr<SDCSolver> solver,
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
        scalar rho,
        scalar delta,
        scalar tol,
        int minSweeps,
        int maxSweeps
        )
        :
        SDC( solver, adaptiveTimeStepper, "gauss-radau", 2, tol, minSweeps, maxSweeps ),
        rho( rho ),
        delta( delta )
    {
        assert( not adaptiveTimeStepper->isEnabled() );
        assert( tol >= delta );
        assert( rho > 0 );

        computeCoefficients();

        k = nodes.rows();

        Info << "Picard Integral Exponential Solver: number of nodes = " << k << endl;

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );

        solver->setNumberOfImplicitStages( k - 1 );
    }

    PIES::PIES(
        scalar rho,
        scalar delta,
        scalar tol
        )
        :
        SDC( "gauss-radau", 2, tol ),
        delta( delta )
    {
        assert( tol >= delta );
        assert( rho > 0 );

        computeCoefficients();

        k = nodes.rows();

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );
    }

    PIES::~PIES()
    {}

    void PIES::computeCoefficients()
    {
        // Initialize variables
        longDouble tstart, tend, tinterval, rho;
        int N, M;

        const longDouble pi = boost::math::constants::pi<longDouble>();

        tstart = 0;
        tend = 1;
        tinterval = tend - tstart;

        // Fine discretization of the semi-disk
        N = 500;
        M = N;

        // Radius of the complex semi-disk S
        rho = longDouble( this->rho );

        // Skeletonization of a semi-disk in the complex plane
        // Discretize the semi-disk with N-steps using polar
        // coordinates. Thereafter, convert the polar coordinates
        // to complex numbers.

        // Compute the number of points on the arc
        longDouble length = longDouble( 0.5 ) * pi * rho + rho;
        longDouble ratio = longDouble( 0.5 ) * pi * rho / length;
        longDouble N_disk_float = ratio * N;

        int N_disk = N_disk_float.convert_to<int>();
        int N_imag = N - N_disk + 1;

        longDouble r = rho;
        longDouble dphi = longDouble( 0.5 ) / (N_disk - 1);
        vector phi( N_disk );
        phi.setZero();

        for ( int i = 0; i < N_disk; i++ )
            phi( i ) = longDouble( 0.5 ) + dphi * i;

        phi *= pi;

        // Convert the polar coordinates r and phi to complex numbers
        vectorc gamma( N );
        gamma.setZero();

        for ( int i = 0; i < phi.rows(); i++ )
            gamma( i ) = std::polar( r, phi( i ) );

        phi.resize( N_imag - 1 );
        phi.setZero();

        longDouble drho = rho / (N_imag - 1);

        for ( int i = 0; i < phi.rows(); i++ )
            phi( i ) = drho * i;

        int index = 0;

        for ( int i = N_disk; i < N_disk + N_imag - 1; i++ )
        {
            gamma( i ) = std::complex<longDouble>( 0, phi( index ) );
            index++;
        }

        // Discretize the time interval
        longDouble dt = tinterval / (M - 1);
        vector t( M );
        t.setZero();

        for ( int i = 0; i < t.rows(); i++ )
            t( i ) = tstart + dt * i;

        // Compute the matrix A
        matrixc A( M, N );

        for ( int i = 0; i < A.rows(); i++ )
            for ( int j = 0; j < A.cols(); j++ )
                A( i, j ) = sdc::exp( gamma( j ) * t( i ) );

        // Matrix compression by QR decomposition with full pivoting
        Eigen::FullPivHouseholderQR<matrixc> qr = A.fullPivHouseholderQr();

        matrixc R = qr.matrixQR().triangularView<Eigen::Upper>();

        int k = 0;

        for ( int i = 0; i < R.rows() - 1; i++ )
        {
            if ( std::abs( R( i + 1, i + 1 ) ) < delta )
            {
                k = i + 1;
                break;
            }
        }

        assert( k > 0 );

        int nbNodes = 0;

        Eigen::VectorXi idx = qr.colsPermutation().indices();

        index = 0;

        for ( int i = 0; i < k; i++ )
        {
            if ( idx( i ) == N_disk or idx( i ) == N_disk - 1 )
                nbNodes += 1;
            else
                nbNodes += 2;
        }

        assert( nbNodes == k * 2 - 2 );

        vectorc gamma_k( nbNodes );

        index = 0;
        int zeroIndex = -1;

        for ( int i = 0; i < k; i++ )
        {
            std::complex<longDouble> value = gamma( idx( i ) );

            if ( idx( i ) == N_disk )
                zeroIndex = index;

            if ( idx( i ) == N_disk || idx( i ) == N_disk - 1 )
            {
                gamma_k( index ) = value;
                index += 1;
            }
            else
            {
                gamma_k( index ) = value;
                gamma_k( index + 1 ) = std::complex<longDouble>( std::real( value ), -std::imag( value ) );
                index += 2;
            }
        }

        assert( zeroIndex > -1 );

        // Compute the coefficients omega

        M = nbNodes;
        dt = tinterval / (M - 1);
        t.resize( M );
        t.setZero();

        for ( int i = 0; i < t.rows(); i++ )
            t( i ) = tstart + dt * i;

        A.resize( nbNodes + 1, M );
        A.setZero();

        for ( int i = 0; i < nbNodes; i++ )
            for ( int j = 0; j < M; j++ )
                A( i, j ) = sdc::exp( gamma_k( i ) * t( j ) );

        // Enforce right-hand-rule based nodes, which results in
        // a L-stable time integration scheme
        A( nbNodes, 0 ) = 1;

        matrixc b( nbNodes + 1, M );
        b.setZero();

        for ( int i = 0; i < nbNodes; i++ )
        {
            for ( int j = 1; j < M; j++ )
            {
                std::complex<longDouble> gamma = gamma_k( i );
                longDouble t0 = t( j - 1 );
                longDouble t1 = t( j );

                if ( i == zeroIndex )
                    b( i, j ) = t1 - t0;
                else
                    b( i, j ) = -( sdc::exp( gamma * t0 ) - sdc::exp( gamma * t1 ) ) / gamma;
            }
        }

        matrixc smatOmega = A.fullPivHouseholderQr().solve( b );
        matrix smatReal = smatOmega.real().transpose();
        matrix smatWeights = smatReal.bottomLeftCorner( smatReal.rows() - 1, smatReal.cols() );
        smatWeights /= tinterval;

        assert( (A * smatOmega - b).norm() < 1.0e-1 * delta );

        // Compute the qmat matrix ( t = 0 .. 1 )

        b.setZero();

        for ( int i = 0; i < nbNodes; i++ )
        {
            for ( int j = 1; j < M; j++ )
            {
                std::complex<longDouble> gamma = gamma_k( i );
                longDouble t0 = tstart;
                longDouble t1 = t( j );

                if ( i == zeroIndex )
                    b( i, j ) = t1 - t0;
                else
                    b( i, j ) = -( sdc::exp( gamma * t0 ) - sdc::exp( gamma * t1 ) ) / gamma;
            }
        }

        matrixc qmatOmega = A.fullPivHouseholderQr().solve( b );
        matrix qmatReal = qmatOmega.real().transpose();
        matrix qmatWeights = qmatReal.bottomLeftCorner( qmatReal.rows() - 1, qmatReal.cols() );
        qmatWeights /= tinterval;

        assert( (A * qmatOmega - b).norm() < 1.0e-1 * delta );

        t.array() += -tstart;
        t /= tinterval;

        nodes = t.cast<scalar>();
        smat = smatWeights.cast<scalar>();
        qmat = qmatWeights.cast<scalar>();
    }
}
