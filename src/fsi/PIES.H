
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *
 * Class
 *  Picard Integral Exponential Solver (PIES)
 */

#ifndef PIES_H
#define PIES_H

#include "SDC.H"
#include <complex>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/constants/constants.hpp>

namespace sdc
{
    namespace mp = boost::multiprecision;

    typedef mp::cpp_dec_float_50 longDouble;
    typedef Eigen::Matrix<std::complex<longDouble>, Eigen::Dynamic, Eigen::Dynamic> matrix50c;
    typedef Eigen::Matrix<longDouble, Eigen::Dynamic, Eigen::Dynamic> matrix50;
    typedef Eigen::Matrix<std::complex<longDouble>, Eigen::Dynamic, 1> vector50c;
    typedef Eigen::Matrix<longDouble, Eigen::Dynamic, 1> vector50;

    template<typename T>
    inline std::complex<T>
    exp( const std::complex<T> & x )
    {
        using std::exp; // use ADL
        T const & r = exp( x.real() );
        return std::polar( r, x.imag() );
    }

    class PIES : public SDC
    {
public:

        PIES(
            std::shared_ptr<SDCSolver> solver,
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
            scalar rho,
            scalar delta,
            scalar tol,
            int minSweeps,
            int maxSweeps
            );

        PIES(
            scalar rho,
            scalar delta,
            scalar tol
            );

        virtual ~PIES();

private:

        void computeCoefficients();

        scalar rho;
        scalar delta;
    };
}


#endif
