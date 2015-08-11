
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
        scalar tol,
        scalar delta
        )
        :
        SDC( solver, adaptiveTimeStepper, "gauss-radau", 2, tol ),
        delta( delta )
    {
        computeCoefficients();

        k = nodes.rows();

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );

        solver->setNumberOfImplicitStages( k - 1 );
    }

    PIES::~PIES()
    {}

    void PIES::computeCoefficients()
    {
        // Initialize variables
        longDouble tstart, tend, tinterval, rho;
        label N, M;

        tstart = -1;
        tend = 1;
        tinterval = tend - tstart;

        // Fine discretization of the semi-disk
        N = 500;

        // Radius of the complex semi-disk S
        rho = 3.15;

        // Skeletonization of a semi-disk in the complex plane
        // Discretize the semi-disk with N-steps using polar
        // coordinates. Thereafter, convert the polar coordinates
        // to complex numbers.

        // Compute the number of points on the arc
        longDouble length = 0.5 * M_PI * rho + rho;
        longDouble ratio = 0.5 * M_PI * rho / length;
        int N_disk = int( ratio * N );
        int N_imag = N - N_disk + 1;

        longDouble r = rho;
        longDouble dphi = 0.5 / (N_disk - 1);
        vector phi( 2 * N_disk - 1 );
        phi.setZero();

        for ( int i = 0; i < N_disk; i++ )
            phi( i ) = 0.5 + dphi * i;

        for ( int i = N_disk; i < phi.rows(); i++ )
            phi( i ) = -phi( i - N_disk );

        phi *= M_PI;

        // Convert the polar coordinates r and phi to complex numbers
        // Gamma exludes the corners of the compex semi-disk
        vectorc gamma( phi.rows() + 2 * N_imag - 1 );
        gamma.setZero();

        for ( int i = 0; i < phi.rows(); i++ )
            gamma( i ) = std::polar( r, phi( i ) );

        phi.resize( 2 * N_imag - 1 );
        phi.setZero();

        for ( int i = 0; i < N_imag; i++ )
            phi( i ) = rho / (N_imag - 1) * i;

        for ( int i = N_imag; i < phi.rows(); i++ )
            phi( i ) = -phi( i - N_imag + 1 );

        int index = 0;

        for ( int i = 2 * N_disk - 1; i < gamma.rows(); i++ )
        {
            gamma( i ) = std::complex<longDouble>( 0, phi( index ) );
            index++;
        }

        N = gamma.rows();
        M = N;

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
                A( i, j ) = std::exp( gamma( j ) * t( i ) );

        // Matrix compression by QR decomposition with full pivoting
        Eigen::FullPivHouseholderQR<matrixc> qr = A.fullPivHouseholderQr();

        matrixc R = qr.matrixQR().triangularView<Eigen::Upper>();

        int nbNodes = 0;

        for ( int i = 0; i < R.rows() - 1; i++ )
        {
            if ( std::abs( R( i + 1, i + 1 ) ) < delta )
            {
                nbNodes = i + 1;
                break;
            }
        }

        vectorc gamma_k( nbNodes );

        Eigen::MatrixXi indices = qr.colsPermutation().indices();

        for ( int i = 0; i < gamma_k.rows(); i++ )
            gamma_k( i ) = gamma( qr.colsPermutation().indices() ( i ) );

        // Compute the coefficients omega

        M = nbNodes;
        dt = tinterval / (M - 1);
        t.resize( M );
        t.setZero();

        for ( int i = 0; i < t.rows(); i++ )
            t( i ) = tstart + dt * i;

        A.resize( nbNodes, M );

        for ( int i = 0; i < A.rows(); i++ )
            for ( int j = 0; j < A.cols(); j++ )
                A( i, j ) = std::abs( gamma_k( i ) * t( j ) );

        matrixc b( A.rows(), A.cols() );
        b.setZero();

        for ( int i = 0; i < b.rows(); i++ )
        {
            for ( int j = 1; j < b.cols(); j++ )
            {
                std::complex<longDouble> gamma = gamma_k( i );
                longDouble t0 = t( j - 1 );
                longDouble t1 = t( j );

                if ( std::abs( gamma ) < 1.0e-14 )
                    b( i, j ) = t1 - t0;
                else
                    b( i, j ) = -( std::exp( gamma * t0 ) - std::exp( gamma * t1 ) ) / gamma;
            }
        }

        matrixc omega = A.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeThinV ).solve( b );
        matrix omegaReal = omega.real();

        // Compute the qmat matrix ( t = -1 .. 1 )

        b.setZero();

        for ( int i = 0; i < b.rows(); i++ )
        {
            for ( int j = 1; j < b.cols(); j++ )
            {
                std::complex<longDouble> gamma = gamma_k( i );
                longDouble t0 = 0;
                longDouble t1 = t( j );

                if ( std::abs( gamma ) < 1.0e-14 )
                    b( i, j ) = t1 - t0;
                else
                    b( i, j ) = -( std::exp( gamma * t0 ) - std::exp( gamma * t1 ) ) / gamma;
            }
        }

        matrixc qmatOmega = A.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeThinV ).solve( b );
        matrix qmatReal = qmatOmega.real();

        t.array() += 1;
        t *= 0.5;
        nodes = t.cast<scalar>();

        smat = omegaReal.cast<scalar>();
        qmat = qmatReal.cast<scalar>();
    }
}
