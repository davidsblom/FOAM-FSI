
#pragma once

/*
 * Take from PFASST library: https://github.com/Parallel-in-Time/PFASST
 */

#include <vector>
using namespace std;


namespace fsi
{
    namespace quadrature
    {
        /**
         * Representation of a polynomial including differentiation, integration and root finding.
         *
         * Nothing more than a \\( n \\)-th order polynomial of the form
         * \\( P_n(x) = \\sum_{i=0}^{n} c_i x^i \\) with coefficients \\( c_i \\),
         * \\( i=0, \\dots, n \\).
         *
         * @tparam CoeffT numerical precision of polynomial coefficients (e.g. `double`)
         */
        template<typename CoeffT>
        class Polynomial
        {
protected:

            /**
             * Coefficients of the polynomial.
             *
             * The coefficient for the highest degree of \\( x \\) has index `0`.
             * The last coefficient is for \\( x^0 \\).
             */
            std::vector<CoeffT> c;

public:

            // ! @{
            Polynomial( size_t n );

            // ! @}

            // ! @{

            /**
             * Order of this polynomial.
             *
             * The order of this polynomial is one less the number of coefficients defined.
             *
             * @returns order of the polynomial
             */
            size_t order() const;

            /**
             * Access coefficient @p i
             *
             * @param[in] i coefficient index (zero-based)
             * @returns \\( i+1 \\)-th coefficient
             * @throws std::out_of_range if index is out of bounds, i.e. @p i >= Polynomial::order()
             *
             * @note The coefficients are stored in reversed order to the degree of the indeterminate,
             *   i.e. `i=0` corresponds to \\( c_n \\) while `i=n` corresponds to \\( c_0 \\).
             *   See also Polynomial::c.
             */
            CoeffT & operator[]( const size_t i );

            /**
             * Differentiate this polynomial.
             *
             * Computes standard differential of this polynomial.
             *
             * @returns differentiated polynomial
             */
            Polynomial<CoeffT> differentiate() const;

            /**
             * Integrates this polynomial.
             *
             * Computes integral of this polynomial.
             *
             * @returns integrated polynomial
             */
            Polynomial<CoeffT> integrate() const;

            // ! @}

            // ! @{

            /**
             * Evaluate polynomial for given value.
             *
             * @tparam xtype numerical type of the value
             * @param[in] x value to evaluate polynomial at
             * @returns value of polynomial at @p x
             */
            template<typename xtype>
            xtype evaluate( const xtype x ) const
            {
                size_t n = this->order();
                xtype v = c[n];

                for ( int j = n - 1; j >= 0; j-- )
                {
                    v = x * v + c[j];
                }

                return v;
            }

            /**
             * Normalizes this polynomial with respect to \\( c_0 \\).
             *
             * @returns normalized polynomial
             */
            Polynomial<CoeffT> normalize() const;

            /**
             * Computes the roots of this polynomial.
             *
             * @returns roots sorted with respect to their value
             */
            std::vector<CoeffT> roots(
                size_t num_iterations = 100,
                CoeffT ztol = 1.0e-20
                ) const;

            // ! @}

            // ! @{

            /**
             * Computes the Legendre polynomial of given order.
             *
             * @param[in] order desired order of the Legendre polynomial
             * @returns Legendre polynomial of order @p order
             */
            static Polynomial<CoeffT> legendre( const size_t order );

            // ! @}
        };
    }
}

#include "Polynomial.tpp"
