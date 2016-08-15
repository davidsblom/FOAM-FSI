
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <cassert>
#include "PetscMatrix.H"
#include "fvCFD.H"

namespace fsi
{
    PetscMatrix::PetscMatrix(
        int rows,
        int cols
        )
        :
        matrix_( new Mat() ),
        rows_( rows ),
        cols_( cols )
    {
        PetscBool petscIsInitialized;
        PetscInitialized( &petscIsInitialized );

        if ( not petscIsInitialized )
        {
            int argc = 1;
            char * arg = new char[9];
            strcpy( arg, "foam-fsi" );
            char ** argv = &arg;

            PetscErrorCode ierr;
            ierr = PetscInitialize( &argc, &argv, "", nullptr );
            CHKERRV( ierr );

            delete [] arg;

            ierr = PetscPushErrorHandler( &PetscMPIAbortErrorHandler, nullptr );
            CHKERRV( ierr );
        }

        PetscErrorCode ierr = 0;
        ierr = MatCreate( PETSC_COMM_WORLD, &*matrix_ );
        CHKERRV( ierr );

        ierr = MatSetType( *matrix_, MATMPIAIJ );
        CHKERRV( ierr );

        ierr = MatSetSizes( *matrix_, PETSC_DECIDE, PETSC_DECIDE, rows, cols );
        CHKERRV( ierr );

        ierr = MatSetUp( *matrix_ );
        CHKERRV( ierr );
    }

    PetscMatrix::PetscMatrix(
        int rows,
        int cols,
        std::unique_ptr<Mat> & matrix
        )
        :
        matrix_( matrix.release() ),
        rows_( rows ),
        cols_( cols )
    {}

    PetscMatrix::PetscMatrix( const PetscMatrix & matrix )
        :
        matrix_( new Mat() ),
        rows_( matrix.rows_ ),
        cols_( matrix.cols_ )
    {
        std::cout << "matrix copy constructor" << std::endl;
    }

    PetscMatrix::~PetscMatrix()
    {
        PetscErrorCode ierr = 0;
        PetscBool petscIsInitialized;
        PetscInitialized( &petscIsInitialized );

        if ( petscIsInitialized )
        {
            // If PetscFinalize is called before ~Matrix
            ierr = MatDestroy( &*matrix_ );
            CHKERRV( ierr );
        }
    }

    void PetscMatrix::compress()
    {
        PetscErrorCode ierr = 0;
        ierr = MatAssemblyBegin( *matrix_, MAT_FINAL_ASSEMBLY );
        CHKERRV( ierr );
        ierr = MatAssemblyEnd( *matrix_, MAT_FINAL_ASSEMBLY );
        CHKERRV( ierr );
    }

    void PetscMatrix::print()
    {
        PetscErrorCode ierr = 0;
        ierr = MatView( *matrix_, PETSC_VIEWER_STDOUT_WORLD );
        CHKERRV( ierr );
    }

    void PetscMatrix::set(
        const int row,
        const int col,
        const PetscScalar value
        )
    {
        assert( row < rows_ );
        assert( col < cols_ );

        PetscErrorCode ierr = 0;
        ierr = MatSetValue( *matrix_, row, col, value, INSERT_VALUES );
        CHKERRV( ierr );
    }

    PetscScalar PetscMatrix::operator()(
        const PetscInt row,
        const PetscInt col
        )
    {
        int size;
        MPI_Comm_size( PETSC_COMM_WORLD, &size );

        if ( size != 1 )
            Foam::abort( FatalError );

        PetscScalar value;

        const PetscErrorCode ierr = MatGetValues( *matrix_, 1, &row, 1, &col, &value );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        return value;
    }

    PetscMatrix operator *(
        const PetscMatrix & A,
        const PetscMatrix & B
        )
    {
        assert( A.cols_ == B.rows_ );

        std::unique_ptr<Mat> C( new Mat() );

        PetscErrorCode ierr = MatMatMult( *(A.matrix_), *(B.matrix_), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &*C );

        PetscMatrix result( A.cols_, B.cols_, C );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        return result;
    }

    PetscVector operator *(
        const PetscMatrix & matrix,
        const PetscVector & vector
        )
    {
        PetscVector result( vector );

        PetscErrorCode ierr = 0;
        ierr = MatMult( *matrix.matrix_, *vector.vector_, *result.vector_ );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        return result;
    }
}
