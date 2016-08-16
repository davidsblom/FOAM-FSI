
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <cassert>
#include "PetscVector.H"
#include "fvCFD.H"

namespace fsi
{
    PetscVector::PetscVector( int rows )
        :
        PetscVector( rows, true )
    {}

    PetscVector::PetscVector(
        int rows,
        bool global
        )
        :
        vector_( new Vec() )
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

        if ( global )
            ierr = VecSetSizes( *vector_, PETSC_DECIDE, rows );
        else
            ierr = VecSetSizes( *vector_, rows, PETSC_DECIDE );

        CHKERRV( ierr );

        ierr = VecSetUp( *vector_ );
        CHKERRV( ierr );
    }

    PetscVector::PetscVector( const PetscVector & vec )
        :
        vector_( new Vec() )
    {
        PetscErrorCode ierr = 0;
        ierr = VecDuplicate( *(vec.vector_), &*vector_ );
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

    void PetscVector::compress()
    {
        PetscErrorCode ierr = 0;
        ierr = VecAssemblyBegin( *vector_ );
        CHKERRV( ierr );
        ierr = VecAssemblyEnd( *vector_ );
        CHKERRV( ierr );
    }

    void PetscVector::print()
    {
        PetscErrorCode ierr = 0;
        ierr = VecView( *vector_, PETSC_VIEWER_STDOUT_WORLD );
        CHKERRV( ierr );
    }

    PetscInt PetscVector::rows()
    {
        PetscInt size;
        PetscErrorCode ierr = 0;
        ierr = VecGetSize( *vector_, &size );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        return size;
    }

    void PetscVector::set(
        const int row,
        const PetscScalar value
        )
    {
        PetscErrorCode ierr = 0;
        ierr = VecSetValue( *vector_, row, value, INSERT_VALUES );
        CHKERRV( ierr );
    }

    PetscScalar PetscVector::operator[]( const PetscInt row )
    {
        int size;
        MPI_Comm_size( PETSC_COMM_WORLD, &size );

        if ( size != 1 )
            Foam::abort( FatalError );

        PetscScalar value;

        PetscErrorCode ierr = VecGetValues( *vector_, 1, &row, &value );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        return value;
    }

    PetscScalar PetscVector::operator()( const PetscInt row )
    {
        return this->operator[]( row );
    }
}
