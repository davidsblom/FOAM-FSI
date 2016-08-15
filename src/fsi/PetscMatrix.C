
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PetscMatrix.H"
#include <cassert>

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

    const Mat & PetscMatrix::get() const
    {
        return *matrix_;
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

    PetscVector operator *(
        const PetscMatrix & matrix,
        const PetscVector & vector
        )
    {
        PetscVector result( vector );

        PetscErrorCode ierr = 0;
        ierr = MatMult( matrix.get(), vector.get(), result.get_mut() );
        assert( ierr == 0 );

        return result;
    }
}
