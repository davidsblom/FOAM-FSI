
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <Eigen/Dense>

namespace fsi
{
    class PetRBFInterpolation
    {
        public:
            PetRBFInterpolation(
                const Eigen::MatrixXd & positions,
                const Eigen::MatrixXd & positionsInterpolation
                );

            ~PetRBFInterpolation();
    };
}
