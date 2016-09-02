
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElRBFInterpolation.H"
#include <cassert>

namespace rbf
{
    ElRBFInterpolation::ElRBFInterpolation()
        :
        H( new El::DistMatrix<double>() ),
        Phi( new El::DistMatrix<double>() )
    {}

    ElRBFInterpolation::ElRBFInterpolation(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        std::unique_ptr<ElDistVector> positions,
        std::unique_ptr<ElDistVector> positionsInterpolation
        )
        :
        H( new El::DistMatrix<double>() ),
        Phi( new El::DistMatrix<double>() )
    {
        compute( std::move( rbfFunction ), std::move( positions ), std::move( positionsInterpolation ) );
    }

    ElRBFInterpolation::~ElRBFInterpolation()
    {}

    void ElRBFInterpolation::compute(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        std::unique_ptr<ElDistVector> positions,
        std::unique_ptr<ElDistVector> positionsInterpolation
        )
    {
        assert( Phi->Height() == 0 );
        assert( Phi->Width() == 0 );
        assert( H->Height() == 0 );
        assert( H->Width() == 0 );

        assert( positions->Width() == positionsInterpolation->Width() );

        H = std::unique_ptr<El::DistMatrix<double> >( new El::DistMatrix<double>( positions->Grid() ) );
        H->AlignWith( *positions );
        El::Zeros( *H, positions->Height(), positions->Height() );

        Phi = std::unique_ptr<El::DistMatrix<double> >( new El::DistMatrix<double>( positionsInterpolation->Grid() ) );
        Phi->AlignRowsWith( *positionsInterpolation );
        El::Zeros( *Phi, positionsInterpolation->Height(), positions->Height() );

        const int dim = positions->Width();

        // Evaluate H matrix

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

        H->Reserve( H->LocalHeight() * H->LocalWidth() );

        for ( int i = 0; i < H->LocalHeight(); i++ )
        {
            const int globalRow = H->GlobalRow( i );

            for ( int j = 0; j < H->LocalWidth(); j++ )
            {
                const int globalCol = H->GlobalCol( j );

                if ( globalCol > globalRow )
                    continue;

                double norm = 0;

                for ( int iDim = 0; iDim < dim; iDim++ )
                    norm += std::pow( pullBuf[i * dim + iDim] - pullBuf[j * dim + iDim + H->LocalHeight() * dim], 2 );

                norm = std::sqrt( norm );

                H->QueueUpdate( globalRow, globalCol, rbfFunction->evaluate( norm ) );
            }
        }

        H->ProcessQueues();

        // Evaluate Phi matrix

        positionsInterpolation->ProcessQueues();

        positions->ReservePulls( Phi->LocalWidth() * dim );
        positionsInterpolation->ReservePulls( Phi->LocalHeight() * dim );

        for ( int i = 0; i < Phi->LocalWidth(); i++ )
            for ( int iDim = 0; iDim < dim; iDim++ )
                positions->QueuePull( Phi->GlobalCol( i ), iDim );

        for ( int i = 0; i < Phi->LocalHeight(); i++ )
            for ( int iDim = 0; iDim < dim; iDim++ )
                positionsInterpolation->QueuePull( Phi->GlobalRow( i ), iDim );

        pullBuf.clear();
        std::vector<double> pullBufInterp;
        positions->ProcessPullQueue( pullBuf );
        positionsInterpolation->ProcessPullQueue( pullBufInterp );

        Phi->Reserve( Phi->LocalHeight() * Phi->LocalWidth() );

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

        Phi->ProcessQueues();
    }

    bool ElRBFInterpolation::initialized()
    {
        return Phi->Height() > 0;
    }

    std::unique_ptr<ElDistVector> ElRBFInterpolation::interpolate( const std::unique_ptr<ElDistVector> & values )
    {
        assert( Phi->Height() > 0 );
        assert( H->Height() > 0 );
        assert( values->Height() == Phi->Width() );

        std::unique_ptr<ElDistVector> result( new ElDistVector( Phi->Grid() ) );

        result->AlignRowsWith( *Phi );

        values->ProcessQueues();

        El::DistMatrix<double> B = *values;

        H->ProcessQueues();

        if ( HCopy.Width() == 0 )
        {
            HCopy = El::DistMatrix<double> ( *H );
            El::DistMatrixReadProxy<double, double, El::MC, El::MR> AProx( HCopy );
            auto & A = AProx.Get();

            p = El::DistPermutation( A.Grid() );
            dSub = El::DistMatrix<double, El::MD, El::STAR> ( A.Grid() );
            El::LDL( A, dSub, p, false, El::LDLPivotCtrl<El::Base<double> >() );
        }

        El::DistMatrixReadProxy<double, double, El::MC, El::MR> AProx( HCopy );
        El::DistMatrixReadWriteProxy<double, double, El::MC, El::MR> BProx( B );
        auto & A = AProx.Get();
        auto & B_LU = BProx.Get();

        El::ldl::SolveAfter( A, dSub, p, B_LU, false );

        Phi->ProcessQueues();

        El::Gemm( El::Orientation::NORMAL, El::Orientation::NORMAL, 1.0, *Phi, B, *result );

        return result;
    }
}
