
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

TEST( Vector, set )
{
    fsi::PetscVector vector( Pstream::nProcs() + 1 );
    vector.set( Pstream::myProcNo(), 2.3 );
    vector.compress();
}

TEST( Vector, copy_constructor )
{
    fsi::PetscVector vector( Pstream::nProcs() + 1 );
    vector.set( Pstream::myProcNo(), 2.3 );
    vector.compress();

    fsi::PetscVector vector2( vector );
    vector2.set( Pstream::myProcNo(), 23.3 );
    vector2.compress();
}
