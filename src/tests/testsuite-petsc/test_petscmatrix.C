
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
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscMatrix matrix( size + 1, size + 1 );
    matrix.set( rank, rank, 5.2 );
    matrix.compress();
}

TEST( Matrix, matmult )
{
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscMatrix matrix( size + 1, size + 1 );
    matrix.set( rank, rank, 5.2 );
    matrix.compress();

    fsi::PetscVector vector1( size + 1 );
    vector1.set( rank, 214.1 );
    vector1.compress();

    fsi::PetscVector result = matrix * vector1;
}
