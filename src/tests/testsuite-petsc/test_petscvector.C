
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
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscVector vector( size + 1 );
    vector.set( rank, 2.3 );
    vector.compress();
}

TEST( Vector, copy_constructor )
{
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscVector vector( size + 1 );
    vector.set( rank, rank );
    vector.compress();

    fsi::PetscVector vector2( vector );
    vector2.set( Pstream::myProcNo(), 23.3 );
    vector2.compress();
}

TEST( Vector, print )
{
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscVector vector( size + 1 );
    vector.set( rank, rank );
    vector.compress();
    vector.print();
}
