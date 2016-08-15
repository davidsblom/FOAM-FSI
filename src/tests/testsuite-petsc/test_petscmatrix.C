
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

TEST( Matrix, view )
{
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscMatrix matrix( size + 1, size + 1 );
    matrix.set( rank, rank, 5.2 );
    matrix.compress();
    matrix.print();
}

TEST( Matrix, matmatmult )
{
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscMatrix A( size + 1, size + 1 );
    A.set( rank, rank, 5.2 );
    A.compress();

    fsi::PetscMatrix B( size + 1, size + 1 );
    B.set( rank, rank, 5.2 );
    B.compress();

    fsi::PetscMatrix C = A * B;

    C.print();
}

TEST( Matrix, get_value )
{
    int rank, size;
    MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    MPI_Comm_size( PETSC_COMM_WORLD, &size );

    fsi::PetscMatrix matrix( size + 1, size + 1 );
    matrix.set( rank, rank, 5.2 );
    matrix.compress();

    if ( size == 1 )
        ASSERT_EQ( matrix( rank, rank ), 5.2 );
}
