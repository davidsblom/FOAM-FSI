
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "PetscMatrix.H"
#include "PetscVector.H"
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

TEST( Matrix, matmult )
{
    fsi::PetscMatrix matrix( Pstream::nProcs() + 1, Pstream::nProcs() + 1 );
    matrix.set( Pstream::myProcNo(), Pstream::myProcNo(), 5.2 );
    matrix.compress();

    fsi::PetscVector vector1( Pstream::nProcs() + 1 );
    vector1.set( Pstream::myProcNo(), 214.1 );
    vector1.compress();

    fsi::PetscVector result = matrix * vector1;
}
