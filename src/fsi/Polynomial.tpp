
/*
 * Copyright [2016] <David Blom>
 */

#include "Polynomial.H"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

namespace fsi
{
    namespace quadrature
    {
        /**
         * @todo Consider issuing a warning/assertion when `n` is zero.
         */
        template<typename CoeffT>
        Polynomial<CoeffT>::Polynomial( std::size_t n )
            : c( n, CoeffT( 0.0 ) )
        {}

        template<typename CoeffT>
        std::size_t Polynomial<CoeffT>::order() const
        {
            return c.size() - 1;
        }

        template<typename CoeffT>
        CoeffT & Polynomial<CoeffT>::operator[]( const std::size_t i )
        {
            return c.at( i );
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::differentiate() const
        {
            Polynomial<CoeffT> p( c.size() - 1 );

            for ( std::size_t j = 1; j < c.size(); j++ )
            {
                p[j - 1] = j * c[j];
            }

            return p;
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::integrate() const
        {
            Polynomial<CoeffT> p( c.size() + 1 );

            for ( std::size_t j = 0; j < c.size(); j++ )
            {
                p[j + 1] = c[j] / (j + 1);
            }

            return p;
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::normalize() const
        {
            Polynomial<CoeffT> p( c.size() );

            for ( std::size_t j = 0; j < c.size(); j++ )
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
            std::size_t num_iterations,
            CoeffT ztol
            ) const
        {
            assert( c.size() >= 1 );
            std::size_t n = c.size() - 1;

            // initial guess
            std::vector<std::complex<CoeffT> > z0( n );

            for ( std::size_t j = 0; j < n; j++ )
            {
                z0[j] = pow( std::complex<CoeffT>( 0.4L, 0.9L ), j );
            }

            // durand-kerner-weierstrass iterations
            Polynomial<CoeffT> p = this->normalize();

            for ( std::size_t k = 0; k < num_iterations; k++ )
            {
                std::complex<CoeffT> num, den;

                for ( std::size_t i = 0; i < n; i++ )
                {
                    num = p.evaluate( z0[i] );
                    den = 1.0L;

                    for ( std::size_t j = 0; j < n; j++ )
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

            for ( std::size_t j = 0; j < n; j++ )
            {
                roots[j] = abs( z0[j] ) < ztol ? 0.0L : real( z0[j] );
            }

            sort( roots.begin(), roots.end() );
            return roots;
        }

        template<typename CoeffT>
        Polynomial<CoeffT> Polynomial<CoeffT>::legendre( const std::size_t order )
        {
            if ( order == 0 )
            {
                Polynomial<CoeffT> p( 1 );
                p[0] = 1.0L;
                return p;
            }

            if ( order == 1 )
            {
                Polynomial<CoeffT> p( 2 );
                p[0] = 0.0L;
                p[1] = 1.0L;
                return p;
            }

            Polynomial<CoeffT> p0( order + 1 ), p1( order + 1 ), p2( order + 1 );
            p0[0] = 1.0L;
            p1[1] = 1.0L;

            // (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
            for ( std::size_t m = 1; m < order; m++ )
            {
                for ( std::size_t j = 1; j < order + 1; j++ )
                {
                    p2[j] = ( (2L * m + 1L) * p1[j - 1] - m * p0[j] ) / (m + 1L);
                }

                p2[0] = -int(m) * p0[0] / (m + 1);

                for ( std::size_t j = 0; j < order + 1; j++ )
                {
                    p0[j] = p1[j];
                    p1[j] = p2[j];
                }
            }

            return p2;
        }

} // namespace quadrature

} //namespace fsi
