
/*
 * Copyright [2016] <David Blom>
 *
 * Class
 *  Picard Integral Exponential Solver (PIES)
 */

#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/constants/constants.hpp>
#include <complex>
#include <numeric>
#include <memory>
#include "SDC.H"

namespace sdc {
template<typename T>
inline std::complex<T>
exp(const std::complex<T> & x) {
    using std::exp;     // use ADL
    T const & r = exp(x.real());
    return std::polar(r, x.imag());
}

template<typename precision>
class PIES : public SDC {
    public:
        PIES(std::shared_ptr<SDCSolver> solver,
            scalar rho,
            scalar delta,
            scalar tol,
            int minSweeps,
            int maxSweeps
            );

        PIES(scalar rho,
            scalar delta,
            scalar tol
            );

        virtual ~PIES();

    private:
        void computeCoefficients();

        scalar rho;
        scalar delta;
};
} // namespace sdc

#include "PIES.tpp"
