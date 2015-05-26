
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include <string>
#include <memory>
#include <ctime>
#include <Eigen/Dense>
#include "RBFFunctionInterface.H"
#include "TPSFunction.H"
#include "fvCFD.H"
#include "debugSwitch.H"

namespace rbf
{
    typedef Eigen::MatrixXd matrix;
    typedef Eigen::VectorXd vector;

    class RBFInterpolation
    {
public:

        RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction );

        void compute(
            const matrix & positions,
            const matrix & positionsInterpolation
            );

        void interpolate(
            const matrix & values,
            matrix & valuesInterpolation
            );

        void interpolate(
            bool polynomialTerm,
            const matrix & positions,
            const matrix & positionsInterpolation,
            const matrix & values,
            matrix & valuesInterpolation
            );

        void interpolate2(
            const matrix & values,
            matrix & valuesInterpolation
            );

        void buildPhi(
            bool polynomialTerm,
            const matrix & positions,
            const matrix & positionsInterpolation
            );

        std::shared_ptr<RBFFunctionInterface> rbfFunction;
        bool computed;
        int n_A;
        int n_B;
        int dimGrid;
        matrix Hhat;
        matrix Phi;
        Eigen::PartialPivLU<matrix> lu;
        Eigen::FullPivLU<matrix> fullPivLu;
        Eigen::LLT<matrix> llt;

        static debug::debugSwitch debug;

private:

        void evaluateH(
            const matrix & positions,
            matrix & H
            );

        void evaluatePhi(
            const matrix & positions,
            const matrix & positionsInterpolation,
            matrix & Phi
            );
    };
}

#endif
