
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

            void set(
                const int row,
                const PetscScalar value
                );

            void view();

            friend PetscVector operator *(
                const PetscMatrix & matrix,
                const PetscVector & vector
                );

        private:
            const std::unique_ptr<Vec> vector_;
            const int rows_;
    };
}
