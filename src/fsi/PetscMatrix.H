
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <petscmat.h>
#include <memory>
#include <iostream>
#include "PetscVector.H"

namespace fsi
{
    class PetscMatrix
    {
        public:
            PetscMatrix(
                int rows,
                int cols
                );

            PetscMatrix(
                int rows,
                int cols,
                std::unique_ptr<Mat> & matrix
                );

            PetscMatrix( const PetscMatrix & matrix );

            ~PetscMatrix();

            void compress();

            void print();

            void set(
                const int row,
                const int col,
                const PetscScalar value
                );

            PetscScalar operator()(
                const PetscInt row,
                const PetscInt col
                );

            friend PetscMatrix operator *(
                const PetscMatrix & A,
                const PetscMatrix & B
                );

            friend PetscVector operator *(
                const PetscMatrix & matrix,
                const PetscVector & vector
                );

        private:
            const std::unique_ptr<Mat> matrix_;

            const int rows_;
            const int cols_;
    };

    PetscMatrix operator *(
        const PetscMatrix & A,
        const PetscMatrix & B
        );

    PetscVector operator *(
        const PetscMatrix & matrix,
        const PetscVector & vector
        );
}
