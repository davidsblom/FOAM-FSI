
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AdaptiveCoarsening.H"

namespace rbf
{
    AdaptiveCoarsening::AdaptiveCoarsening(
        double tol,
        double reselectionTol,
        int minPoints,
        int maxPoints
        ) :
        tol( tol ),
        reselectionTol( reselectionTol ),
        minPoints( minPoints ),
        maxPoints( maxPoints ),
        rbf( new ElRBFInterpolation() )
    {
        assert( maxPoints >= minPoints );
        assert( tol > 0 );
        assert( tol <= 1 );
    }

    AdaptiveCoarsening::~AdaptiveCoarsening()
    {}

    void AdaptiveCoarsening::compute(
        std::shared_ptr<RBFFunctionInterface> function,
        std::unique_ptr<El::DistMatrix<double> > pos,
        std::unique_ptr<El::DistMatrix<double> > posInterpolation
        )
    {
        // Selection is performed as necessary. It is only performed when interpolate()
        // is called.
        pos->ProcessQueues();
        posInterpolation->ProcessQueues();

        // Store the untouched data for later reselection of points
        this->rbfFunction = function;
        this->positions = std::move( pos );
        this->positionsInterpolation = std::move( posInterpolation );
    }

    std::pair<int, double> AdaptiveCoarsening::computeError( const std::unique_ptr<El::DistMatrix<double> > & values )
    {
        // Select a subset of values based on the selected points

        std::unique_ptr<El::DistMatrix<double> > valuesCoarse( new El::DistMatrix<double>() );
        El::Zeros( *valuesCoarse, selectedPositions.size(), values->Width() );
        selectData( values, valuesCoarse );

        // Perform the interpolation

        std::unique_ptr<El::DistMatrix<double> > result = rbfCoarse->interpolate( valuesCoarse );

        assert( values->Height() == result->Height() );

        // Compute error
        El::DistMatrix<double> diff = *values;
        El::Axpy( -1, *result, diff );
        El::DistMatrix<double, El::MC, El::STAR> errors;
        El::RowTwoNorms( diff, errors );

        // Get location of max error
        El::Entry<double> locMax = El::MaxAbsLoc( errors );

        // Scale by largest value
        El::RowTwoNorms( *values, errors );
        double maxValue = El::MaxAbs( errors );

        if ( maxValue != 0 )
            return std::pair<int, double>( locMax.i, locMax.value / maxValue );
        else
            return std::pair<int, double>( locMax.i, locMax.value );
    }

    void AdaptiveCoarsening::greedySelection( const std::unique_ptr<El::DistMatrix<double> > & values )
    {
        selectedPositions.clear();

        selectedPositions.push_back( 0 );
        selectedPositions.push_back( 1 );
        int maxPoints = std::min( this->maxPoints, positions->Height() );
        int minPoints = std::min( this->minPoints, positions->Height() );
        double error = 0;

        for ( int i = 0; i < maxPoints; i++ )
        {
            // Build the matrices for the RBF interpolation
            std::unique_ptr<El::DistMatrix<double> > positionsCoarse( new El::DistMatrix<double>() );
            El::Zeros( *positionsCoarse, selectedPositions.size(), positions->Width() );
            selectData( positions, positionsCoarse );

            // Perform the RBF interpolation
            std::unique_ptr<El::DistMatrix<double> > positionsInterpolationCoarse( new El::DistMatrix<double>() );
            *positionsInterpolationCoarse = *positions;
            rbfCoarse = std::unique_ptr<ElRBFInterpolation>( new ElRBFInterpolation( rbfFunction, std::move( positionsCoarse ), std::move( positionsInterpolationCoarse ) ) );

            std::pair<int, double> largestError = computeError( values );
            error = largestError.second;

            // Break if maximum points are reached
            if ( int( selectedPositions.size() ) >= maxPoints )
                break;

            bool convergence = error < tol && int( selectedPositions.size() ) >= minPoints;

            if ( convergence )
                break;

            selectedPositions.push_back( largestError.first );
        }

        if ( El::mpi::Rank() == 0 )
        {
            std::cout << "RBF interpolation coarsening: selected ";
            std::cout << selectedPositions.size() << "/" << positions->Height() << " points,";
            std::cout << " error = " << error << ", tol = " << tol << std::endl;
        }

        // Initialize interpolator

        std::unique_ptr<El::DistMatrix<double> > positionsCoarse( new El::DistMatrix<double>() );
        El::Zeros( *positionsCoarse, selectedPositions.size(), positions->Width() );

        selectData( positions, positionsCoarse );

        rbf = std::unique_ptr<ElRBFInterpolation>( new ElRBFInterpolation() );

        std::unique_ptr<El::DistMatrix<double> > positionsInterpolationTmp( new El::DistMatrix<double>( *positionsInterpolation ) );

        rbf->compute( rbfFunction, std::move( positionsCoarse ), std::move( positionsInterpolationTmp ) );
    }

    bool AdaptiveCoarsening::initialized()
    {
        return rbf->initialized();
    }

    std::unique_ptr<El::DistMatrix<double> > AdaptiveCoarsening::interpolate( const std::unique_ptr<El::DistMatrix<double> > & values )
    {
        values->ProcessQueues();

        bool greedyPerformed = false;

        // Greedy selection never performed => do it now
        if ( not rbf->initialized() )
        {
            // Only do a greedy selection if the values actually mean something
            double maxAbs = El::MaxAbs( *values );

            if ( maxAbs > 0 )
            {
                greedySelection( values );
                greedyPerformed = true;
            }
            else
            {
                std::unique_ptr<El::DistMatrix<double> > result( new El::DistMatrix<double>() );
                El::Zeros( *result, positionsInterpolation->Height(), positionsInterpolation->Width() );
                return result;
            }
        }

        // Evaluate error
        std::pair<int, double> largestError = computeError( values );

        if ( !greedyPerformed && El::mpi::Rank() == 0 )
        {
            std::cout << "RBF interpolation coarsening: ";
            std::cout << " error = " << largestError.second << ", tol = " << reselectionTol;
            std::cout << ", reselection = ";

            if ( largestError.second >= reselectionTol )
                std::cout << "true";
            else
                std::cout << "false";

            std::cout << std::endl;
        }

        if ( largestError.second >= reselectionTol && !greedyPerformed )
        {
            // the error is too large. Do a reselection
            greedySelection( values );
        }

        std::unique_ptr<El::DistMatrix<double> > selectedValues( new El::DistMatrix<double>() );
        El::Zeros( *selectedValues, selectedPositions.size(), values->Width() );

        selectData( values, selectedValues );

        return rbf->interpolate( selectedValues );
    }

    void AdaptiveCoarsening::selectData(
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
