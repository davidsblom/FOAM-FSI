
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PetscVector.H"
#include <cassert>

namespace fsi
{
    PetscVector::PetscVector( int rows )
        :
        vector_( new Vec() ),
        rows_( rows )
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
        ierr = VecCreate( PETSC_COMM_WORLD, &*vector_ );
        CHKERRV( ierr );

        ierr = VecSetType( *vector_, VECMPI );
        CHKERRV( ierr );

        ierr = VecSetSizes( *vector_, PETSC_DECIDE, rows );
        CHKERRV( ierr );

        ierr = VecSetUp( *vector_ );
        CHKERRV( ierr );
    }

    PetscVector::~PetscVector()
    {
        PetscErrorCode ierr = 0;
        PetscBool petscIsInitialized;
        PetscInitialized( &petscIsInitialized );

        if ( petscIsInitialized )
        {
            // If PetscFinalize is called before ~Matrix
            ierr = VecDestroy( &*vector_ );
            CHKERRV( ierr );
        }
    }
}
