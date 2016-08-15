
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
    class PetscVector
    {
        public:
            PetscVector( int rows );

            PetscVector( const PetscVector & vec );

            ~PetscVector();

            const Vec & get() const;

            Vec & get_mut();

            void compress();

            void set(
                const int row,
                const PetscScalar value
                );

        private:
            const std::unique_ptr<Vec> vector_;
            const int rows_;
    };
}
