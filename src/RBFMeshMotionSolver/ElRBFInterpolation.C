
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElRBFInterpolation.H"
#include <cassert>

namespace rbf
{
    ElRBFInterpolation::ElRBFInterpolation(
        std::unique_ptr<RBFFunctionInterface> function,
        std::unique_ptr<El::DistMultiVec<double> > pos,
        std::unique_ptr<El::DistMultiVec<double> > posInterp
        )
        :
        rbfFunction( std::move( function ) ),
        positions( std::move( pos ) ),
        positionsInterpolation( std::move( posInterp ) ),
        H( new El::DistMatrix<double>() )
    {
        assert( positions->Width() == positionsInterpolation->Width() );

        El::Zeros( *H, positions->Height(), positions->Height() );

        const int dim = positions->Width();

        // Evaluate H matrix

        H->Reserve( H->LocalHeight() * H->LocalWidth() );

        positions->ProcessQueues();
        positionsInterpolation->ProcessQueues();

        for ( int i = 0; i < H->LocalHeight(); i++ )
        {
            const int globalRow = H->GlobalRow( i );

            for ( int j = 0; j < H->LocalWidth(); j++ )
            {
                double norm = 0;
                const int globalCol = H->GlobalCol( j );

                for ( int iDim = 0; iDim < dim; iDim++ )
                    norm += std::pow( positions->Get( globalRow, iDim ) - positions->Get( globalCol, iDim ), 2 );

                norm = std::sqrt( norm );

                H->QueueUpdate( globalRow, globalCol, rbfFunction->evaluate( norm ) );
            }
        }
    }

    ElRBFInterpolation::~ElRBFInterpolation()
    {}

    std::unique_ptr<El::DistMatrix<double> > ElRBFInterpolation::interpolate( const std::unique_ptr<El::DistMatrix<double> > values )
    {
        std::unique_ptr<El::DistMatrix<double> > result( new El::DistMatrix<double>() );
        const int dim = positions->Width();

        // Evaluate Phi matrix

        El::DistMatrix<double> Phi( positionsInterpolation->Height(), positions->Height() );

        Phi.Reserve( Phi.LocalHeight() * Phi.LocalWidth() );

        for ( int i = 0; i < Phi.LocalHeight(); i++ )
        {
            const int globalRow = Phi.GlobalRow( i );

            for ( int j = 0; j < Phi.LocalWidth(); j++ )
            {
                double norm = 0;
                const int globalCol = Phi.GlobalCol( j );

                for ( int iDim = 0; iDim < dim; iDim++ )
                    norm += std::pow( positions->Get( globalCol, iDim ) - positionsInterpolation->Get( globalRow, iDim ), 2 );

                norm = std::sqrt( norm );

                Phi.QueueUpdate( globalRow, globalCol, rbfFunction->evaluate( norm ) );
            }
        }

        values->ProcessQueues();

        El::DistMatrix<double> B = *values;

        H->ProcessQueues();

        El::LinearSolve( *H, B );

        Phi.ProcessQueues();

        El::Gemm( El::Orientation::NORMAL, El::Orientation::NORMAL, 1.0, Phi, B, *result );

        return result;
    }
}
