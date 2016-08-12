
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PetscMatrix.H"

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
        }

        PetscErrorCode ierr = 0;
        ierr = MatCreate( PETSC_COMM_WORLD, &*matrix_ );
        CHKERRV( ierr );

        ierr = MatSetType( *matrix_, MATSEQAIJ );
        CHKERRV( ierr );

        ierr = MatSetSizes( *matrix_, PETSC_DECIDE, PETSC_DECIDE, rows, cols );
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
}
