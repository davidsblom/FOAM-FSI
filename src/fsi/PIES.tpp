
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *
 * This class implements the Picard Integral Exponential Solver as proposed
 * in D. Kushnir and V. Rokhlin, "A Highly Accurate Solver for Stiff Ordinary
 * Differential Equations" (2011).
 * The basis functionality is implemented in the SDC class. The method
 * updateCoefficients() computes the coefficients for the PIES solver.
 */

#include <math.h>
#include "GaussRadau.H"

namespace sdc
{
    template<typename precision>
    PIES<precision>::PIES(
        std::shared_ptr<SDCSolver> solver,
        scalar rho,
        scalar delta,
        scalar tol,
        int minSweeps,
        int maxSweeps
        )
        :
        SDC( solver, std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>(2) ), tol, minSweeps, maxSweeps ),
        rho( rho ),
        delta( delta )
    {
        assert( rho > 0 );

        std::clock_t start;
        double duration;

        start = std::clock();

        computeCoefficients();

        duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

        k = nodes.rows();

        Info << "Picard Integral Exponential Solver: number of nodes = " << k << ", timing = " << duration << " s" << Foam::endl;

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );

        solver->setNumberOfImplicitStages( k - 1 );

        data->initialize( k, solver->getDOF() );
    }

    template<typename precision>
    PIES<precision>::PIES(
        scalar rho,
        scalar delta,
        scalar tol
        )
        :
        SDC( std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>(2) ), tol ),
        rho( rho ),
        delta( delta )
    {
        assert( rho > 0 );

        computeCoefficients();

        k = nodes.rows();

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );
    }

    template<typename precision>
    PIES<precision>::~PIES()
    {}

    template<typename precision>
    void PIES<precision>::computeCoefficients()
    {
        typedef Eigen::Matrix<std::complex<precision>, Eigen::Dynamic, Eigen::Dynamic> matrix50c;
        typedef Eigen::Matrix<precision, Eigen::Dynamic, Eigen::Dynamic> matrix50;
        typedef Eigen::Matrix<std::complex<precision>, Eigen::Dynamic, 1> vector50c;
        typedef Eigen::Matrix<precision, Eigen::Dynamic, 1> vector50;

        Info << "PIES: compute integration coefficients" << Foam::endl;

        // Initialize variables
        precision tstart, tend, tinterval, rho, delta, eps;
        int N, M;

        const precision pi = boost::math::constants::pi<precision>();

        tstart = -1;
        tend = 1;
        tinterval = tend - tstart;

        // Fine discretization of the semi-disk
        N = 800;
        M = N;

        // Radius of the complex semi-disk S
        rho = precision( this->rho );

        // Selection criteria of nodes on the semi-desik S
        eps = 1.0e-18;
        delta = precision( this->delta );

        // Skeletonization of a semi-disk in the complex plane

        // Compute the number of points on the arc
        precision length = precision( 0.5 ) * pi * rho + rho;
        precision dx = length / N;

        // Travel from the origin over the semi circle

        vector50c gamma( N + 3 );
        gamma.setZero();
        gamma( N ) = std::complex<precision>( -rho, 0 );
        gamma( N + 1 ) = std::complex<precision>( 0, -rho );
        gamma( N + 2 ) = std::complex<precision>( 0, rho );

        Info << "PIES: discretize gamma on semi-disk" << Foam::endl;

        for ( int i = 0; i < N; i++ )
        {
            precision distanceTravelled = i * dx;

            if ( distanceTravelled <= rho )
                gamma( i ) = std::complex<precision>( 0, distanceTravelled );
            else
            {
                using std::cos;
                using std::sin;
                precision distanceTravelledOnDisk = distanceTravelled - rho + 0.5 * M_PI * rho;
                precision theta = distanceTravelledOnDisk / rho;
                precision x = cos( theta ) * rho;
                precision y = sin( theta ) * rho;
                gamma( i ) = std::complex<precision>( x, y );
            }
        }

        // Discretize the time interval
        precision dt = tinterval / (M - 1);
        vector50 t( M );
        t.setZero();

        for ( int i = 0; i < t.rows(); i++ )
            t( i ) = tstart + dt * i;

        Info << "PIES: build matrix A of size " << t.rows() << " x " << gamma.rows() << Foam::endl;

        // Compute the matrix A
        matrix50c A( t.rows(), gamma.rows() );

        for ( int i = 0; i < A.rows(); i++ )
            for ( int j = 0; j < A.cols(); j++ )
                A( i, j ) = sdc::exp( gamma( j ) * t( i ) );

        Info << "PIES: compute QR decomposition of A" << Foam::endl;

        // Matrix compression by QR decomposition with full pivoting
        Eigen::FullPivHouseholderQR<matrix50c> qr = A.fullPivHouseholderQr();

        matrix50c R = qr.matrixQR().template triangularView<Eigen::Upper>();

        int k = 0;

        for ( int i = 0; i < R.rows() - 1; i++ )
        {
            k = i;

            if ( std::abs( R( i + 1, i + 1 ) ) / std::abs( R( 0, 0 ) ) < delta )
            {
                k = i + 1;
                break;
            }
        }

        assert( k > 0 );

        int nbNodes = 0;

        Eigen::VectorXi idx = qr.colsPermutation().indices();

        for ( int i = 0; i < k; i++ )
        {
            if ( idx( i ) == 0 or idx( i ) == N or idx( i ) == N + 1 or idx( i ) == N + 2 )
                nbNodes += 1;
            else
                nbNodes += 2;
        }

        vector50c gamma_k( nbNodes );

        int index = 0;
        int zeroIndex = -1;

        Info << "PIES: build subset lambda" << Foam::endl;

        for ( int i = 0; i < k; i++ )
        {
            std::complex<precision> value = gamma( idx( i ) );

            if ( idx( i ) == 0 )
                zeroIndex = index;

            if ( idx( i ) == 0 or idx( i ) == N or idx( i ) == N + 1 or idx( i ) == N + 2 )
            {
                gamma_k( index ) = value;
                index += 1;
            }
            else
            {
                gamma_k( index ) = value;
                gamma_k( index + 1 ) = std::complex<precision>( std::real( value ), -std::imag( value ) );
                index += 2;
            }
        }

        // Compute the coefficients omega

        // One extra node due to adding of constraint
        nbNodes += 1;

        Info << "PIES: compute the coeffients omega" << Foam::endl;

        for ( int iter = 0; iter < 20; iter++ )
        {
            M = nbNodes;
            dt = tinterval / (M - 1);
            t.resize( M );
            t.setZero();

            for ( int i = 0; i < t.rows(); i++ )
                t( i ) = tstart + dt * i;

            A.resize( gamma_k.rows() + 1, M );
            A.setZero();

            // Ensure that the system is never over-determined
            assert( A.rows() <= A.cols() );

            for ( int i = 0; i < gamma_k.rows(); i++ )
                for ( int j = 0; j < M; j++ )
                    A( i, j ) = sdc::exp( gamma_k( i ) * t( j ) );

            // Enforce right-hand-rule based nodes, which results in
            // a L-stable time integration scheme
            A( gamma_k.rows(), 0 ) = 1;

            matrix50c b( A.rows(), A.cols() );
            b.setZero();

            for ( int i = 0; i < gamma_k.rows(); i++ )
            {
                for ( int j = 1; j < M; j++ )
                {
                    std::complex<precision> gamma = gamma_k( i );
                    precision t0 = t( j - 1 );
                    precision t1 = t( j );

                    if ( i == zeroIndex )
                        b( i, j ) = t1 - t0;
                    else
                        b( i, j ) = -( sdc::exp( gamma * t0 ) - sdc::exp( gamma * t1 ) ) / gamma;
                }
            }

            matrix50c smatOmega = A.fullPivHouseholderQr().solve( b );
            matrix50 smatReal = smatOmega.real().transpose();
            matrix50 smatWeights = smatReal.bottomLeftCorner( smatReal.rows() - 1, smatReal.cols() );
            smatWeights /= tinterval;

            // Compute the qmat matrix ( t = 0 .. 1 )

            b.setZero();

            for ( int i = 0; i < gamma_k.rows(); i++ )
            {
                for ( int j = 1; j < M; j++ )
                {
                    std::complex<precision> gamma = gamma_k( i );
                    precision t0 = tstart;
                    precision t1 = t( j );

                    if ( i == zeroIndex )
                        b( i, j ) = t1 - t0;
                    else
                        b( i, j ) = -( sdc::exp( gamma * t0 ) - sdc::exp( gamma * t1 ) ) / gamma;
                }
            }

            matrix50c qmatOmega = A.fullPivHouseholderQr().solve( b );
            matrix50 qmatReal = qmatOmega.real().transpose();
            matrix50 qmatWeights = qmatReal.bottomLeftCorner( qmatReal.rows() - 1, qmatReal.cols() );
            qmatWeights /= tinterval;

            t.array() += -tstart;
            t /= tinterval;

            nodes = t.template cast<scalar>();
            smat = smatWeights.template cast<scalar>();
            qmat = qmatWeights.template cast<scalar>();

            precision residual = (A * qmatOmega - b).norm();
            std::cout << "PIES: residual = " << residual;
            std::cout << ", number of nodes = " << nbNodes;
            std::cout << ", delta = " << delta;
            std::cout << std::endl;

            if ( residual < eps )
                break;

            nbNodes++;
        }
    }
}
