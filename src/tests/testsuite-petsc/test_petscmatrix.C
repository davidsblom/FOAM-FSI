
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "PetscMatrix.H"

TEST( Matrix, constructor )
{
    fsi::PetscMatrix matrix( 5, 5 );
}
