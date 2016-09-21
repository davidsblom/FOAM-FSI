
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "ElRBFInterpolation.H"
#include "TPSFunction.H"
#include "RBFFunctionInterface.H"
#include "NoCoarsening.H"

using namespace rbf;

TEST( NoCoarsener, interpolate )
{
    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positions( new El::DistMatrix<double, El::VR, El::STAR>() );
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positionsInterpolation( new El::DistMatrix<double, El::VR, El::STAR>() );
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > data( new El::DistMatrix<double, El::VR, El::STAR>() );
    El::Zeros( *positionsInterpolation, 8, 1 );
    El::Zeros( *positions, 4, 1 );
    El::Zeros( *data, positions->Height(), 1 );

    const double dx = 1.0 / positions->Height();
    const double dy = 1.0 / double( positionsInterpolation->Height() );

    positions->Reserve( positions->LocalHeight() * positions->LocalWidth() );

    for ( int i = 0; i < positions->LocalHeight(); i++ )
    {
        const int globalRow = positions->GlobalRow( i );

        for ( int j = 0; j < positions->LocalWidth(); j++ )
        {
            const int globalCol = positions->GlobalCol( j );
            positions->QueueUpdate( globalRow, globalCol, (globalRow + 1) * dx );
        }
    }

    positionsInterpolation->Reserve( positionsInterpolation->LocalHeight() * positionsInterpolation->LocalWidth() );

    for ( int i = 0; i < positionsInterpolation->LocalHeight(); i++ )
    {
        const int globalRow = positionsInterpolation->GlobalRow( i );

        for ( int j = 0; j < positionsInterpolation->LocalWidth(); j++ )
        {
            const int globalCol = positionsInterpolation->GlobalCol( j );
            positionsInterpolation->QueueUpdate( globalRow, globalCol, (globalRow + 1) * dy );
        }
    }

    ElRBFInterpolation rbf( rbfFunction, std::move( positions ), std::move( positionsInterpolation ) );

    // Interpolate some data

    data->Reserve( data->LocalHeight() * data->LocalWidth() );

    for ( int i = 0; i < data->LocalHeight(); i++ )
    {
        const int globalRow = data->GlobalRow( i );

        for ( int j = 0; j < data->LocalWidth(); j++ )
        {
            const int globalCol = data->GlobalCol( j );
            const double x = dx * (globalCol + 1);

            data->QueueUpdate( globalRow, globalCol, std::sin( x ) );
        }
    }

    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > result = rbf.interpolate( std::move( data ) );
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > result2;

    {
        std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positions( new El::DistMatrix<double, El::VR, El::STAR>() );
        std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positionsInterpolation( new El::DistMatrix<double, El::VR, El::STAR>() );
        std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > data( new El::DistMatrix<double, El::VR, El::STAR>() );
        El::Zeros( *positionsInterpolation, 8, 1 );
        El::Zeros( *positions, 4, 1 );
        El::Zeros( *data, positions->Height(), 1 );

        const double dx = 1.0 / positions->Height();
        const double dy = 1.0 / double( positionsInterpolation->Height() );

        positions->Reserve( positions->LocalHeight() * positions->LocalWidth() );

        for ( int i = 0; i < positions->LocalHeight(); i++ )
        {
            const int globalRow = positions->GlobalRow( i );

            for ( int j = 0; j < positions->LocalWidth(); j++ )
            {
                const int globalCol = positions->GlobalCol( j );
                positions->QueueUpdate( globalRow, globalCol, (globalRow + 1) * dx );
            }
        }

        positionsInterpolation->Reserve( positionsInterpolation->LocalHeight() * positionsInterpolation->LocalWidth() );

        for ( int i = 0; i < positionsInterpolation->LocalHeight(); i++ )
        {
            const int globalRow = positionsInterpolation->GlobalRow( i );

            for ( int j = 0; j < positionsInterpolation->LocalWidth(); j++ )
            {
                const int globalCol = positionsInterpolation->GlobalCol( j );
                positionsInterpolation->QueueUpdate( globalRow, globalCol, (globalRow + 1) * dy );
            }
        }

        NoCoarsening rbf;
        rbf.compute( rbfFunction, std::move( positions ), std::move( positionsInterpolation ) );

        // Interpolate some data

        data->Reserve( data->LocalHeight() * data->LocalWidth() );

        for ( int i = 0; i < data->LocalHeight(); i++ )
        {
            const int globalRow = data->GlobalRow( i );

            for ( int j = 0; j < data->LocalWidth(); j++ )
            {
                const int globalCol = data->GlobalCol( j );
                const double x = dx * (globalCol + 1);

                data->QueueUpdate( globalRow, globalCol, std::sin( x ) );
            }
        }

        result2 = rbf.interpolate( std::move( data ) );
    }

    std::vector<double> buffer1, buffer2;
    result->ReservePulls( result->Height() * result->Width() );
    result2->ReservePulls( result2->Height() * result2->Width() );

    for ( int i = 0; i < result->Height(); i++ )
        for ( int j = 0; j < result->Width(); j++ )
            result->QueuePull( i, j );

    for ( int i = 0; i < result2->Height(); i++ )
        for ( int j = 0; j < result2->Width(); j++ )
            result2->QueuePull( i, j );

    result->ProcessPullQueue( buffer1 );
    result2->ProcessPullQueue( buffer2 );

    EXPECT_EQ( buffer1.size(), buffer2.size() );
    EXPECT_EQ( result->Height(), result2->Height() );
    EXPECT_EQ( result->Width(), result2->Width() );

    for ( size_t i = 0; i < buffer1.size(); i++ )
        EXPECT_NEAR( buffer1[i], buffer2[i], 1e-14 );
}
