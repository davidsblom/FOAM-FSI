
/*
* Author
*   David Blom, TU Delft. All rights reserved.
*/

#ifndef inverseMatrixLibrary_H
#define inverseMatrixLibrary_H

#include <string>
#include <memory>
#include <Eigen/Dense>
#include "RBFFunctionInterface.H"

namespace rbf
{
    typedef Eigen::MatrixXd matrix;
    typedef Eigen::VectorXd vector;

    class inverseMatrixLibrary
    {
    public:

        inverseMatrixLibrary( std::shared_ptr<RBFFunctionInterface> rbfFunction );

        void updateInverseSymmetry(
            const matrix & oldCinverse,
            const matrix & newpoint,
            const matrix & oldPositions,
            matrix & Cinverse
        );

        void updateInverse(
            const matrix & oldCinverse,
            const matrix & newpoint,
            const matrix & oldPositions,
            matrix & Cinverse
        );

        void computeInverseEigen(
            const matrix & positions,
            matrix & Cinverse
        );

        void solveSystemEigen(
            const matrix & positions,
            const matrix & displacement,
            matrix & solution
        );

        /*void compute(
        const matrix & positions,
        const matrix & positionsInterpolation
        );

        void interpolate(
            const matrix & values,
            matrix & valuesInterpolation
        );

        void interpolate(
            const matrix & positions,
            const matrix & positionsInterpolation,
            const matrix & values,
            matrix & valuesInterpolation
        );*/

        void getControlMatrix(
            const matrix & positions,
            matrix & Crbf
        );

        void getControlMatrixLower(
            const matrix & positions,
            matrix & Crbf
        );

        void getControlMatrixUpper(
            const matrix & positions,
            matrix & Crbf
        );


    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    int n_A;
    int dimGrid;
    matrix C;
    };
}

#endif
