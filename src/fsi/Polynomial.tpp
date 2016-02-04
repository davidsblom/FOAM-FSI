
#include "Polynomial.H"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
using namespace std;


namespace fsi
{
    namespace quadrature
    {
        /**
         * @todo Consider issuing a warning/assertion when `n` is zero.
         */
        template<typename CoeffT>
        Polynomial<CoeffT>::Polynomial( size_t n )
            : c( n, CoeffT( 0.0 ) )
        {}

        template<typename CoeffT>
        size_t Polynomial<CoeffT>::order() const
        {
            return c.size() - 1;
        }

        template<typename CoeffT>
        CoeffT & Polynomial<CoeffT>::operator[]( const size_t i )
        {
            return c.at( i );
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::differentiate() const
        {
            Polynomial<CoeffT> p( c.size() - 1 );

            for ( size_t j = 1; j < c.size(); j++ )
            {
                p[j - 1] = j * c[j];
            }

            return p;
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::integrate() const
        {
            Polynomial<CoeffT> p( c.size() + 1 );

            for ( size_t j = 0; j < c.size(); j++ )
            {
                p[j + 1] = c[j] / (j + 1);
            }

            return p;
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::normalize() const
        {
            Polynomial<CoeffT> p( c.size() );

            for ( size_t j = 0; j < c.size(); j++ )
            {
                p[j] = c[j] / c.back();
            }

            return p;
        }

        /**
         * @internals
         * @note Asserts this polynomial has at least order 1 if `NDEBUG` is not defined.
         * @endinternals
         */
        template<typename CoeffT>
        std::vector<CoeffT> Polynomial<CoeffT>::roots(
            size_t num_iterations,
            CoeffT ztol
            ) const
        {
            assert( c.size() >= 1 );
            size_t n = c.size() - 1;

            // initial guess
            std::vector<complex<CoeffT> > z0( n );

            for ( size_t j = 0; j < n; j++ )
            {
                z0[j] = pow( complex<double>( 0.4, 0.9 ), j );
            }

            // durand-kerner-weierstrass iterations
            Polynomial<CoeffT> p = this->normalize();

            for ( size_t k = 0; k < num_iterations; k++ )
            {
                complex<CoeffT> num, den;

                for ( size_t i = 0; i < n; i++ )
                {
                    num = p.evaluate( z0[i] );
                    den = 1.0;

                    for ( size_t j = 0; j < n; j++ )
                    {
                        if ( j == i )
                        {
                            continue;
                        }

                        den = den * (z0[i] - z0[j]);
                    }

                    z0[i] = z0[i] - num / den;
                }
            }

            std::vector<CoeffT> roots( n );

            for ( size_t j = 0; j < n; j++ )
            {
                roots[j] = abs( z0[j] ) < ztol ? 0.0 : real( z0[j] );
            }

            sort( roots.begin(), roots.end() );
            return roots;
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::legendre( const size_t order )
        {
            if ( order == 0 )
            {
                Polynomial<CoeffT> p( 1 );
                p[0] = 1.0;
                return p;
            }

            if ( order == 1 )
            {
                Polynomial<CoeffT> p( 2 );
                p[0] = 0.0;
                p[1] = 1.0;
                return p;
            }

            Polynomial<CoeffT> p0( order + 1 ), p1( order + 1 ), p2( order + 1 );
            p0[0] = 1.0;
            p1[1] = 1.0;

            // (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
            for ( size_t m = 1; m < order; m++ )
            {
                for ( size_t j = 1; j < order + 1; j++ )
                {
                    p2[j] = ( (2 * m + 1) * p1[j - 1] - m * p0[j] ) / (m + 1);
                }

                p2[0] = -int(m) * p0[0] / (m + 1);

                for ( size_t j = 0; j < order + 1; j++ )
                {
                    p0[j] = p1[j];
                    p1[j] = p2[j];
                }
            }

            return p2;
        }
    }
}
