
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
                bool global
                );

            PetscMatrix( std::unique_ptr<Mat> & matrix );

            PetscMatrix( const PetscMatrix & matrix );

            ~PetscMatrix();

            PetscInt cols() const;

            void compress() const;

            void print() const;

            PetscInt rows() const;

            void set(
                const int row,
                const int col,
                const PetscScalar value
                );

            void setLocal(
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

            friend PetscVector solve(
                const PetscMatrix & mat,
                const PetscVector & vec
                );

        private:
            const std::unique_ptr<Mat> matrix_;
    };

    PetscMatrix operator *(
        const PetscMatrix & A,
        const PetscMatrix & B
        );

    PetscVector operator *(
        const PetscMatrix & matrix,
        const PetscVector & vector
        );

    PetscVector solve(
        const PetscMatrix & mat,
        const PetscVector & vec
        );
}
