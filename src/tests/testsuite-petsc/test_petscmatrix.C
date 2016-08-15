
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "PetscMatrix.H"
#include "fvCFD.H"

TEST( Matrix, constructor )
{
    fsi::PetscMatrix matrix( 5, 5 );
}

TEST( Matrix, set )
{
    fsi::PetscMatrix matrix( Pstream::nProcs() + 1, Pstream::nProcs() + 1 );
    matrix.set( Pstream::myProcNo(), Pstream::myProcNo(), 5.2 );
    matrix.compress();
}
