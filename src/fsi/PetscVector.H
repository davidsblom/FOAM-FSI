
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

            PetscVector(
                int rows,
                bool global
                );

            PetscVector( const PetscVector & vec );

            ~PetscVector();

            void compress();

            void print();

            PetscInt rows();

            void set(
                const int row,
                const PetscScalar value
                );

            void setLocal(
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
    };
}
