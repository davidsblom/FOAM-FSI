
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "PetscVector.H"
#include "fvCFD.H"

TEST( Vector, constructor )
{
    fsi::PetscVector vector( 5 );
}
