
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <mxx/sort.hpp>
#include "UnitCoarsening.H"

namespace rbf
{
    UnitCoarsening::UnitCoarsening(
        double tol,
        int minPoints,
        int maxPoints
        )
        :
        tol( tol ),
        minPoints( minPoints ),
        maxPoints( maxPoints )
    {
        assert( maxPoints >= minPoints );
        assert( tol > 0 );
        assert( tol <= 1 );
    }

    UnitCoarsening::~UnitCoarsening()
    {}

    void UnitCoarsening::compute(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        std::unique_ptr<El::DistMatrix<double> > positions,
        std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
        )
    {
        // selection is only performed once. Thereafter, the selected points are used for the
        // interpolation.

        selectedPositions.clear();

        selectedPositions.push_back( 0 );
        selectedPositions.push_back( 1 );
        int maxPoints = std::min( this->maxPoints, positions->Height() );
        int minPoints = std::min( this->minPoints, positions->Height() );

        double largestError = 0;

        for ( int i = 0; i < maxPoints; i++ )
        {
            // Build the matrices for the RBF interpolation
            std::unique_ptr<El::DistMatrix<double> > positionsCoarse( new El::DistMatrix<double>() );
            std::unique_ptr<El::DistMatrix<double> > valuesCoarse( new El::DistMatrix<double>() );
            El::Ones( *valuesCoarse, selectedPositions.size(), positions->Width() );
            El::Zeros( *positionsCoarse, selectedPositions.size(), positions->Width() );

            selectData( positions, positionsCoarse );

            // Perform the RBF interpolation
            std::unique_ptr<El::DistMatrix<double> > positionsInterpolationCoarse( new El::DistMatrix<double>() );

            *positionsInterpolationCoarse = *positions;

            ElRBFInterpolation rbf( rbfFunction, std::move( positionsCoarse ), std::move( positionsInterpolationCoarse ) );

            std::unique_ptr<El::DistMatrix<double> > result = rbf.interpolate( std::move( valuesCoarse ) );

            assert( result->Height() == positions->Height() );
            assert( result->Width() == positions->Width() );

            // Evaluate the error locally
            result->ReservePulls( result->LocalHeight() * result->Width() );

            for ( int k = 0; k < result->LocalHeight(); k++ )
                for ( int l = 0; l < result->Width(); l++ )
                    result->QueuePull( result->GlobalRow( k ), l );

            std::vector<double> buffer;
            result->ProcessPullQueue( buffer );

            std::vector<double> localErrors( result->LocalHeight() );

            int index = 0;

            for ( int k = 0; k < result->LocalHeight(); k++ )
            {
                double norm = 0;

                for ( int l = 0; l < result->Width(); l++ )
                {
                    norm += std::pow( std::abs( buffer[index] - 1.0 ), 2 );
                    index++;
                }

                localErrors[k] = std::sqrt( norm );
            }

            // Locally, select the point with the largest error
            auto it = std::max_element( localErrors.begin(), localErrors.end() );
            int largestErrorLocalIndex = std::distance( localErrors.begin(), it );
            int largestErrorGlobalIndex = result->GlobalRow( largestErrorLocalIndex );
            largestError = -1;

            // Add safeguard in case the number cpus is larger than positions->Height()
            if ( it != localErrors.end() )
                largestError = *it;

            // Send the largest error, and corresponding index to all (!) cpus.
            std::vector<double> globalErrors = mxx::allgather( largestError );
            std::vector<int> largestErrorGlobalIndexes = mxx::allgather( largestErrorGlobalIndex );

            // Select the point with the largest error of the global distribution
            it = std::max_element( globalErrors.begin(), globalErrors.end() );
            int largestErrorIndex = std::distance( globalErrors.begin(), it );
            largestError = *it;
            int globalIndexError = largestErrorGlobalIndexes[largestErrorIndex];

            // Break if maximum points are reached
            if ( int( selectedPositions.size() ) >= maxPoints )
                break;

            bool convergence = largestError < tol && int( selectedPositions.size() ) >= minPoints;

            if ( convergence )
                break;

            selectedPositions.push_back( globalIndexError );
        }

        if ( El::mpi::Rank() == 0 )
        {
            std::cout << "RBF interpolation coarsening: selected ";
            std::cout << selectedPositions.size() << "/" << positions->Height() << " points,";
            std::cout << " error = " << largestError << ", tol = " << tol << std::endl;
        }

        std::unique_ptr<El::DistMatrix<double> > positionsCoarse( new El::DistMatrix<double>() );
        El::Zeros( *positionsCoarse, selectedPositions.size(), positions->Width() );

        selectData( positions, positionsCoarse );

        rbf.compute( rbfFunction, std::move( positionsCoarse ), std::move( positionsInterpolation ) );
    }

    bool UnitCoarsening::initialized()
    {
        return rbf.initialized();
    }

    std::unique_ptr<El::DistMatrix<double> > UnitCoarsening::interpolate( const std::unique_ptr<El::DistMatrix<double> > & values )
    {
        std::unique_ptr<El::DistMatrix<double> > selectedValues( new El::DistMatrix<double>() );
        El::Zeros( *selectedValues, selectedPositions.size(), values->Width() );

        selectData( values, selectedValues );

        return rbf.interpolate( selectedValues );
    }

    void UnitCoarsening::selectData(
        const std::unique_ptr<El::DistMatrix<double> > & data,
        std::unique_ptr<El::DistMatrix<double> > & selection
        )
    {
        assert( selection->Height() == int( selectedPositions.size() ) );

        data->ProcessQueues();

        int nbPulls = 0;

        for ( size_t j = 0; j < selectedPositions.size(); j++ )
        {
            for ( int iDim = 0; iDim < data->Width(); iDim++ )
            {
                if ( selection->Owner( j, iDim ) == El::mpi::Rank( selection->DistComm() ) )
                    nbPulls++;
            }
        }

        data->ReservePulls( nbPulls );

        for ( size_t j = 0; j < selectedPositions.size(); j++ )
        {
            for ( int iDim = 0; iDim < data->Width(); iDim++ )
            {
                if ( selection->Owner( j, iDim ) == El::mpi::Rank( selection->DistComm() ) )
                    data->QueuePull( selectedPositions[j], iDim );
            }
        }

        std::vector<double> buffer;
        data->ProcessPullQueue( buffer );
        int index = 0;

        for ( size_t j = 0; j < selectedPositions.size(); j++ )
        {
            for ( int iDim = 0; iDim < data->Width(); iDim++ )
            {
                if ( selection->Owner( j, iDim ) == El::mpi::Rank( selection->DistComm() ) )
                {
                    selection->Set( j, iDim, buffer[index] );
                    index++;
                }
            }
        }
    }
}
