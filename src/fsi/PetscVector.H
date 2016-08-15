
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <petscvec.h>
#include <memory>
#include <iostream>

namespace fsi
{
    class PetscMatrix;

    class PetscVector
    {
        public:
            PetscVector( int rows );

            PetscVector( const PetscVector & vec );

            ~PetscVector();

            void compress();

            void print();

            void set(
                const int row,
                const PetscScalar value
                );

            PetscScalar operator[]( const PetscInt row );

            PetscScalar operator()( const PetscInt row );

            friend PetscVector operator *(
                const PetscMatrix & matrix,
                const PetscVector & vector
                );

            friend PetscVector solve(
                const PetscMatrix & mat,
                const PetscVector & vec
                );

        private:
            const std::unique_ptr<Vec> vector_;
            const int rows_;
    };
}
