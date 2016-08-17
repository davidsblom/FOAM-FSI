
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <cassert>
#include "PetscMatrix.H"
#include "fvCFD.H"
#include <petscksp.h>

namespace fsi
{
    PetscMatrix::PetscMatrix(
        int rows,
        int cols
        )
        :
        PetscMatrix( rows, cols, true )
    {}

    PetscMatrix::PetscMatrix(
        int rows,
        int cols,
        bool global
        )
        :
        matrix_( new Mat() )
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

        if ( global )
            ierr = MatSetSizes( *matrix_, PETSC_DECIDE, PETSC_DECIDE, rows, cols );
        else
            ierr = MatSetSizes( *matrix_, rows, cols, PETSC_DECIDE, PETSC_DECIDE );

        CHKERRV( ierr );

        ierr = MatSetUp( *matrix_ );
        CHKERRV( ierr );
    }

    PetscMatrix::PetscMatrix( std::unique_ptr<Mat> & matrix )
        :
        matrix_( matrix.release() )
    {}

    PetscMatrix::PetscMatrix( const PetscMatrix & matrix )
        :
        matrix_( new Mat() )
    {
        PetscErrorCode ierr = 0;
        ierr = MatDuplicate( *(matrix.matrix_), MAT_COPY_VALUES, &*matrix_ );
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

    PetscInt PetscMatrix::cols() const
    {
        PetscInt rows, cols;
        PetscErrorCode ierr = 0;
        ierr = MatGetSize( *matrix_, &rows, &cols );

        if ( ierr != 1 )
            Foam::abort( FatalError );

        return cols;
    }

    void PetscMatrix::compress() const
    {
        PetscErrorCode ierr = 0;
        ierr = MatAssemblyBegin( *matrix_, MAT_FINAL_ASSEMBLY );
        CHKERRV( ierr );
        ierr = MatAssemblyEnd( *matrix_, MAT_FINAL_ASSEMBLY );
        CHKERRV( ierr );
    }

    void PetscMatrix::print() const
    {
        PetscErrorCode ierr = 0;
        ierr = MatView( *matrix_, PETSC_VIEWER_STDOUT_WORLD );
        CHKERRV( ierr );
    }

    PetscInt PetscMatrix::rows() const
    {
        PetscInt rows, cols;
        PetscErrorCode ierr = 0;
        ierr = MatGetSize( *matrix_, &rows, &cols );

        if ( ierr != 1 )
            Foam::abort( FatalError );

        return rows;
    }

    void PetscMatrix::set(
        const int row,
        const int col,
        const PetscScalar value
        )
    {
        assert( row < rows() );
        assert( col < cols() );

        PetscErrorCode ierr = 0;
        ierr = MatSetValue( *matrix_, row, col, value, INSERT_VALUES );
        CHKERRV( ierr );
    }

    void PetscMatrix::setLocal(
        const int row,
        const int col,
        const PetscScalar value
        )
    {
        PetscErrorCode ierr = 0;
        ierr = MatSetValuesLocal( *matrix_, 1, &row, 1, &col, &value, INSERT_VALUES );
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
        assert( A.cols() == B.rows() );

        std::unique_ptr<Mat> C( new Mat() );

        PetscErrorCode ierr = MatMatMult( *(A.matrix_), *(B.matrix_), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &*C );

        PetscMatrix result( C );

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

    PetscVector solve(
        const PetscMatrix & matrix,
        const PetscVector & vector
        )
    {
        PetscErrorCode ierr = 0;

        KSP ksp;
        ierr = KSPCreate( PETSC_COMM_WORLD, &ksp );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        ierr = KSPSetOperators( ksp, *(matrix.matrix_), *(matrix.matrix_) );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        ierr = KSPSetType( ksp, KSPGMRES );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        ierr = KSPSetUp( ksp );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        PetscVector result( vector );

        ierr = KSPSolve( ksp, *(vector.vector_), *(result.vector_) );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        KSPConvergedReason reason;

        ierr = KSPGetConvergedReason( ksp, &reason );

        if ( ierr != 0 )
            Foam::abort( FatalError );

        if ( reason == KSP_DIVERGED_INDEFINITE_PC )
        {
            std::cout << "Divergence because of indefinite preconditioner" << std::endl;
            std::cout << "Run the executable again but with -pc_factor_shift_positive_definite option" << std::endl;
            Foam::abort( FatalError );
        }
        else
        if ( reason < 0 )
        {
            std::cout << "Other kind of divergence: this should not happen." << std::endl;
            Foam::abort( FatalError );
        }
        else
        {
            PetscInt its;
            KSPGetIterationNumber( ksp, &its );
            std::cout << "Convergence in " << its << " iterations." << std::endl;
        }

        return result;
    }
}
