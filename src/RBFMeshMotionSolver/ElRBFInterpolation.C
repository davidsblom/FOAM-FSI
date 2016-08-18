
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElRBFInterpolation.H"
#include <cassert>

namespace rbf
{
    ElRBFInterpolation::ElRBFInterpolation(
        std::unique_ptr<RBFFunctionInterface> rbfFunction,
        std::unique_ptr<El::DistMatrix<double> > positions,
        std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
        )
        :
        H( new El::DistMatrix<double>() ),
        Phi( new El::DistMatrix<double>() )
    {
        assert( positions->Width() == positionsInterpolation->Width() );

        El::Zeros( *H, positions->Height(), positions->Height() );
        El::Zeros( *Phi, positionsInterpolation->Height(), positions->Height() );

        const int dim = positions->Width();

        // Evaluate H matrix

        H->Reserve( H->LocalHeight() * H->LocalWidth() );

        positions->ProcessQueues();

        positions->ReservePulls( H->LocalHeight() * positions->Width() + H->LocalWidth() * positions->Width() );

        for ( int i = 0; i < H->LocalHeight(); i++ )
        {
            for ( int iDim = 0; iDim < dim; iDim++ )
                positions->QueuePull( H->GlobalRow( i ), iDim );
        }

        for ( int j = 0; j < H->LocalWidth(); j++ )
        {
            for ( int iDim = 0; iDim < dim; iDim++ )
                positions->QueuePull( H->GlobalCol( j ), iDim );
        }

        std::vector<double> pullBuf;
        positions->ProcessPullQueue( pullBuf );

        for ( int i = 0; i < H->LocalHeight(); i++ )
        {
            const int globalRow = H->GlobalRow( i );

            for ( int j = 0; j < H->LocalWidth(); j++ )
            {
                double norm = 0;
                const int globalCol = H->GlobalCol( j );

                for ( int iDim = 0; iDim < dim; iDim++ )
                    norm += std::pow( pullBuf[i * dim + iDim] - pullBuf[j * dim + iDim + H->LocalHeight() * dim], 2 );

                norm = std::sqrt( norm );

                H->QueueUpdate( globalRow, globalCol, rbfFunction->evaluate( norm ) );
            }
        }

        // Evaluate Phi matrix

        Phi->Reserve( Phi->LocalHeight() * Phi->LocalWidth() );

        positionsInterpolation->ProcessQueues();

        positions->ReservePulls( Phi->LocalWidth() * dim );
        positionsInterpolation->ReservePulls( Phi->LocalHeight() * dim );

        for ( int i = 0; i < Phi->LocalWidth(); i++ )
            for ( int iDim = 0; iDim < dim; iDim++ )
                positions->QueuePull( Phi->GlobalRow( i ), iDim );

        for ( int i = 0; i < Phi->LocalHeight(); i++ )
            for ( int iDim = 0; iDim < dim; iDim++ )
                positionsInterpolation->QueuePull( Phi->GlobalRow( i ), iDim );

        pullBuf.clear();
        std::vector<double> pullBufInterp;
        positions->ProcessPullQueue( pullBuf );
        positionsInterpolation->ProcessPullQueue( pullBufInterp );

        for ( int i = 0; i < Phi->LocalHeight(); i++ )
        {
            const int globalRow = Phi->GlobalRow( i );

            for ( int j = 0; j < Phi->LocalWidth(); j++ )
            {
                double norm = 0;
                const int globalCol = Phi->GlobalCol( j );

                for ( int iDim = 0; iDim < dim; iDim++ )
                    norm += std::pow( pullBuf[j * dim + iDim] - pullBufInterp[i * dim + iDim], 2 );

                norm = std::sqrt( norm );

                Phi->QueueUpdate( globalRow, globalCol, rbfFunction->evaluate( norm ) );
            }
        }
    }

    ElRBFInterpolation::~ElRBFInterpolation()
    {}

    std::unique_ptr<El::DistMatrix<double> > ElRBFInterpolation::interpolate( const std::unique_ptr<El::DistMatrix<double> > values )
    {
        std::unique_ptr<El::DistMatrix<double> > result( new El::DistMatrix<double>() );

        values->ProcessQueues();

        El::DistMatrix<double> B = *values;

        H->ProcessQueues();

        El::LinearSolve( *H, B );

        Phi->ProcessQueues();

        El::Gemm( El::Orientation::NORMAL, El::Orientation::NORMAL, 1.0, *Phi, B, *result );

        return result;
    }
}
