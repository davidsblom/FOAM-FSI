
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

namespace sdc
{
    typedef long double longDouble;
    typedef Eigen::Matrix<std::complex<longDouble>, Eigen::Dynamic, Eigen::Dynamic> matrixc;
    typedef Eigen::Matrix<longDouble, Eigen::Dynamic, Eigen::Dynamic> matrix;
    typedef Eigen::Matrix<std::complex<longDouble>, Eigen::Dynamic, 1> vectorc;
    typedef Eigen::Matrix<longDouble, Eigen::Dynamic, 1> vector;

    class PIES : public SDC
    {
public:

        PIES(
            std::shared_ptr<SDCSolver> solver,
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
            scalar tol,
            scalar delta
            );

        virtual ~PIES();

private:

        void computeCoefficients();

        scalar delta;
    };
}


#endif
