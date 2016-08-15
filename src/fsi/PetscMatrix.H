
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <petscmat.h>
#include <memory>
#include <iostream>

namespace fsi
{
    class PetscMatrix
    {
        public:
            PetscMatrix(
                int rows,
                int cols
                );

            ~PetscMatrix();

            void compress();

            void set(
                const int row,
                const int col,
                const PetscScalar value
                );

        private:
            const std::unique_ptr<Mat> matrix_;

            const int rows_;
            const int cols_;
    };
}
