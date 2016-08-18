
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "ElRBFInterpolation.H"
#include "TPSFunction.H"
#include "RBFInterpolation.H"

using namespace rbf;

TEST( ElRBFInterpolation, constructor )
{
    std::unique_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    std::unique_ptr<El::DistMatrix<double> > positions( new El::DistMatrix<double>( 4, 1 ) );
    std::unique_ptr<El::DistMatrix<double> > positionsInterpolation( new El::DistMatrix<double>( 8, 1 ) );

    const double dx = 1.0 / (positions->Height() - 1);
    const double dy = 1.0 / (positionsInterpolation->Height() - 1);

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

    ElRBFInterpolation rbf( std::move( rbfFunction ), std::move( positions ), std::move( positionsInterpolation ) );
}

TEST( ElRBFInterpolation, interpolate )
{
    std::unique_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    std::unique_ptr<El::DistMatrix<double> > positions( new El::DistMatrix<double>() );
    std::unique_ptr<El::DistMatrix<double> > positionsInterpolation( new El::DistMatrix<double>() );
    std::unique_ptr<El::DistMatrix<double> > data( new El::DistMatrix<double>() );
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

    ElRBFInterpolation rbf( std::move( rbfFunction ), std::move( positions ), std::move( positionsInterpolation ) );

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

    std::unique_ptr<El::DistMatrix<double> > result = rbf.interpolate( std::move( data ) );

    for ( int i = 0; i < result->Height(); i++ )
    {
        const double x = dx * (i + 1);
        const double error = std::abs( std::sin( x ) - result->Get( i, 0 ) );
        EXPECT_LT( error, 1.0 );
    }
}

TEST( ElRBFInterpolation, verify_Eigen )
{
    std::unique_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    std::unique_ptr<El::DistMatrix<double> > positions( new El::DistMatrix<double>() );
    std::unique_ptr<El::DistMatrix<double> > positionsInterpolation( new El::DistMatrix<double>() );
    std::unique_ptr<El::DistMatrix<double> > data( new El::DistMatrix<double>() );
    El::Zeros( *positionsInterpolation, 8, 1 );
    El::Zeros( *positions, 4, 1 );
    El::Zeros( *data, positions->Height(), 1 );
    matrix positionsEigen( positions->Height(), positions->Width() ), positionsInterpolationEigen( positionsInterpolation->Height(), positionsInterpolation->Width() );
    std::shared_ptr<RBFFunctionInterface> rbfFunctionEigen( new TPSFunction() );
    matrix valuesEigen( positionsEigen.rows(), positionsEigen.cols() ), valuesInterpolationEigen;

    const double dx = 1.0 / positions->Height();
    const double dy = 1.0 / double( positionsInterpolation->Height() );

    for ( int i = 0; i < positionsEigen.rows(); i++ )
        positionsEigen( i, 0 ) = (i + 1) * dx;

    for ( int i = 0; i < positionsInterpolationEigen.rows(); i++ )
        positionsInterpolationEigen( i, 0 ) = ( i + 1 ) * dy;

    for ( int i = 0; i < valuesEigen.rows(); i++ )
        valuesEigen( i, 0 ) = std::sin( (i+1) * dx );

    RBFInterpolation rbfEigen( rbfFunctionEigen, false, true );
    rbfEigen.compute( positionsEigen, positionsInterpolationEigen );
    rbfEigen.interpolate( valuesEigen, valuesInterpolationEigen );

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

    positions->ProcessQueues();
    positionsInterpolation->ProcessQueues();

    ElRBFInterpolation rbf( std::move( rbfFunction ), std::move( positions ), std::move( positionsInterpolation ) );

    // Interpolate some data

    data->Reserve( data->LocalHeight() * data->LocalWidth() );

    for ( int i = 0; i < data->LocalHeight(); i++ )
    {
        const int globalRow = data->GlobalRow( i );

        for ( int j = 0; j < data->LocalWidth(); j++ )
        {
            const int globalCol = data->GlobalCol( j );
            const double x = dx * (globalRow + 1);

            data->QueueUpdate( globalRow, globalCol, std::sin( x ) );
        }
    }

    El::Display( *data );
    std::cout << valuesEigen << std::endl;

    std::unique_ptr<El::DistMatrix<double> > result = rbf.interpolate( std::move( data ) );

    El::Display( *result );

    std::cout << valuesInterpolationEigen << std::endl;

}
