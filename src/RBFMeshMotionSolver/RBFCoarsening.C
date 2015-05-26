
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFCoarsening.H"

namespace rbf
{
    Foam::debug::debugSwitch RBFCoarsening::debug( "RBFInterpolation", 0 );

    RBFCoarsening::RBFCoarsening( std::shared_ptr<RBFInterpolation> rbf )
        :
        rbf( rbf ),
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction ) ) ),
        enabled( false ),
        livePointSelection( false ),
        livePointSelectionSumValues( false ),
        tol( 0 ),
        tolLivePointSelection( 0 ),
        coarseningMinPoints( 0 ),
        coarseningMaxPoints( 0 ),
        twoPointSelection ( false ),
        exportTxt( false ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 )
    {
        assert( rbf );
    }

    RBFCoarsening::RBFCoarsening(
        std::shared_ptr<RBFInterpolation> rbf,
        bool enabled,
        bool livePointSelection,
        bool livePointSelectionSumValues,
        double tol,
        double tolLivePointSelection,
        int coarseningMinPoints,
        int coarseningMaxPoints,
        bool twoPointSelection,
        bool exportTxt
        )
        :
        rbf( rbf ),
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction ) ) ),
        enabled( enabled ),
        livePointSelection( livePointSelection ),
        livePointSelectionSumValues( livePointSelectionSumValues ),
        tol( tol ),
        tolLivePointSelection( tolLivePointSelection ),
        coarseningMinPoints( coarseningMinPoints ),
        coarseningMaxPoints( coarseningMaxPoints ),
        twoPointSelection ( twoPointSelection ),
        exportTxt( exportTxt ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 )
    {
        assert( rbf );
        assert( coarseningMinPoints <= coarseningMaxPoints );
        assert( coarseningMinPoints > 0 );
        assert( coarseningMaxPoints > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( tolLivePointSelection > 0 );
        assert( tolLivePointSelection < 1 );
    }

    /* Select a subset of control point with a greedy algorithm.
     * The selection of the points is based on a displacement/motion of
     * a unit displacement at every control point. Based on a user specified tolerance,
     * the control points are selected.
     */
    void RBFCoarsening::compute(
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        this->positions = positions;
        this->positionsInterpolation = positionsInterpolation;
    }

    void RBFCoarsening::greedySelection( const matrix & values )
    {
        assert( positions.cols() == positionsInterpolation.cols() );
        assert( positions.cols() > 0 );
        assert( positions.rows() > 0 );
        assert( positionsInterpolation.rows() > 0 );
        assert( positions.rows() == values.rows() );

        matrix usedPositions = positions;

        if ( enabled )
        {
            // Greedy algorithm

            rbf::vector errorList( positions.rows() );
            selectedPositions.resize( 2 );

            for ( int i = 0; i < selectedPositions.rows(); i++ )
                selectedPositions( i ) = i;

            assert( positions.rows() >= selectedPositions.rows() );

            rbf::matrix positionsInterpolationCoarse = positions;

            int maxNbPoints = std::min( coarseningMaxPoints, static_cast<int>( positions.rows() ) );
            int minPoints = std::min( coarseningMinPoints, static_cast<int>( positions.rows() ) );
            double error = 0;

            // Create RBF interpolator

            // Initialize the matrix Phi if more than 1 positions have been selected
            // before the greedening algorithm starts
            rbfCoarse->Phi.resize( 0, 0 );

            /*Info << "0: selectedPositions.rows() = " << selectedPositions.rows() << endl;
            for ( int i = 1; i < selectedPositions.rows(); i++ )
            {
                rbf::matrix positionsCoarse( i, positions.cols() );
                Info << "0: i = " << i << endl;
                for ( int j = 0; j < i; j++ )
                {
                    Info << "0: j = " << i << endl;
                    positionsCoarse.row( j ) = positions.row( selectedPositions( j ) );
                }

                rbfCoarse->buildPhi( livePointSelection, positionsCoarse, positionsInterpolationCoarse );
            }*/


            // Run the greedy algorithm
            double runTimeInterpolate = 0.0;
            double runTimeError = 0.0;
            double runTimeConvergence = 0.0;
            bool addedSecondPoint = false;
            int counter = selectedPositions.rows();

            Info << "twoPointSelection = " << twoPointSelection << endl;
            while(true)
            {
                std::clock_t t = std::clock();

                // Build the matrices used for the RBF interpolation
                rbf::matrix positionsCoarse( counter, positions.cols() );
                rbf::matrix valuesCoarse( positionsCoarse.rows(), positionsCoarse.cols() );
                rbf::matrix valuesInterpolationCoarse( positionsInterpolationCoarse.rows(), positionsInterpolationCoarse.cols() );

                for ( int j = 0; j < selectedPositions.rows(); j++ )
                {
                    positionsCoarse.row( j ) = positions.row( selectedPositions( j ) );
                    valuesCoarse.row( j ) = values.row( selectedPositions( j ) );
                }

                // Perform the RBF interpolation.
                rbfCoarse->interpolate( livePointSelection, positionsCoarse, positionsInterpolationCoarse, valuesCoarse, valuesInterpolationCoarse );

                if ( debug > 0 )
                {
                    t = std::clock() - t;
                    runTimeInterpolate += static_cast<float>(t) / CLOCKS_PER_SEC;
                    t = std::clock();
                }


                // Evaluate the error
                for ( int j = 0; j < valuesInterpolationCoarse.rows(); j++ )
                    errorList( j ) = ( valuesInterpolationCoarse.row( j ).array() - values.row( j ).array() ).matrix().norm();



                // Select the point with the largest error which is not already selected.
                int index = -1;
                double largestError = -1;
                for ( int j = 0; j < errorList.rows(); j++ )
                {
                    // Check if points are selected or not. Is quiet costly (Ns x Nc) and
                    // should not be needed since error on selected points should be zero (and thus not maximum)
                    // However, for the case where there is zero displacement all errors are zero.
                    // Additional function to check whether the largestError = 0 (<SMALL) and do select next consecutive point

                    /*bool notSelected = true;

                    for ( int k = 0; k < selectedPositions.rows(); k++ )
                        if ( selectedPositions( k ) == j )
                            notSelected = false;

                    if ( errorList( j ) > largestError && notSelected )*/
                    if ( errorList( j ) > largestError )
                    {
                        index = j;
                        largestError = errorList( j );
                    }
                }


                // Additional function to check whether the largestError = 0 (<SMALL) and do select next consecutive point
                if(largestError<SMALL)
                {
                    index = selectedPositions.rows();
                }

                int index2 = -1;
                double largestError2 = -1;
                //selected point with largest error in opposite direction (more than 90 degrees differenc in direction)
                if(twoPointSelection)
                {
                    vector largestErrorVector = valuesInterpolationCoarse.row( index ) - values.row( index );
                    for ( int j = 0; j < errorList.rows(); j ++ )
                    {
                        vector errorVector = valuesInterpolationCoarse.row( j ) - values.row( j );
                        if(largestErrorVector.dot(errorVector)<-SMALL && largestError2 < errorList(j) )
                        {
                            index2 = j;
                            largestError2 = errorList( j );
                        }
                    }
                }

                if ( debug > 0 )
                {
                    t = std::clock() - t;
                    runTimeError += static_cast<float>(t) / CLOCKS_PER_SEC;
                    t = std::clock();
                    //Info << index << ": largestError = " << largestError << ", " << index2 << ": largestError2 = " << largestError2 << endl;
                }

                double epsilon = std::sqrt( SMALL );
                error = ( errorList ).matrix().norm() / (values.norm() + epsilon);
                bool convergence = (error < tol && counter >= minPoints) || counter >= maxNbPoints;

                if ( convergence )
                {
                    break;
                }

                selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
                selectedPositions( selectedPositions.rows() - 1 ) = index;
                counter ++;

                //Break if maximum point are reached
                if ( counter >= maxNbPoints )
                {
                    break;
                }

                //Add second point if possible
                if(twoPointSelection && index2 >= 0 && index != index2 )
                {
                    addedSecondPoint = true;
                    selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
                    selectedPositions( selectedPositions.rows() - 1 ) = index2;
                    counter ++;
                }

                if ( debug > 0 )
                {
                    t = std::clock() - t;
                    runTimeConvergence += static_cast<float>(t) / CLOCKS_PER_SEC;
                    t = std::clock();
                }
            }
            if ( debug > 0 )
            {
                Info << "RBFCoarsening::debug 1. interpolate to surface = " << runTimeInterpolate << " s" << endl;
                Info << "RBFCoarsening::debug 2. find largest error = " << runTimeError << " s" <<". Added second point = " << addedSecondPoint <<  endl;
                Info << "RBFCoarsening::debug 3. convergence check = " << runTimeConvergence << " s" << endl;
            }

            Info << "RBF interpolation coarsening: selected " << selectedPositions.rows() << "/" << positions.rows() << " points, error = " << error << ", tol = " << tol << endl;

            rbf::matrix positionsCoarse( selectedPositions.rows(), positions.cols() );

            for ( int i = 0; i < selectedPositions.rows(); i++ )
                positionsCoarse.row( i ) = positions.row( selectedPositions( i ) );

            usedPositions = positionsCoarse;

            if ( exportTxt )
            {
                std::string fileNameTXT = "rbf-coarsening-greedy-selection-" + std::to_string( fileExportIndex ) + ".txt";
                std::ofstream fileTXT( fileNameTXT );

                if ( fileTXT.is_open() )
                    fileTXT << usedPositions;

                std::string fileNameCSV = "rbf-coarsening-greedy-selection-" + std::to_string( fileExportIndex ) + ".csv";
                std::ofstream fileCSV( fileNameCSV );
                Eigen::IOFormat CSVFmt( Eigen::FullPrecision, Eigen::DontAlignCols, "," );

                if ( fileCSV.is_open() )
                    fileCSV << usedPositions.format( CSVFmt );

                fileExportIndex++;
            }
        }

        rbf->compute( usedPositions, positionsInterpolation );
    }

    void RBFCoarsening::interpolate(
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        matrix usedValues = values;

        if ( enabled )
        {
            if ( livePointSelection )
            {
                // For RBF mesh interpolation, the values to be interpolated need to be
                // the total displacements. As input, the incremental displacements
                // are given.
                if ( livePointSelectionSumValues )
                {
                    if ( this->values.cols() != values.cols() )
                        this->values = values;
                    else
                        this->values.array() += values.array();
                }
                else
                    this->values = values;

                // Check if re-selection is necessary
                bool reselection = true;

                if ( rbfCoarse->computed )
                {
                    rbf::matrix valuesCoarse( selectedPositions.rows(), values.cols() );
                    rbf::matrix valuesInterpolationCoarse( positions.rows(), valuesInterpolation.cols() );

                    for ( int j = 0; j < selectedPositions.rows(); j++ )
                        valuesCoarse.row( j ) = this->values.row( selectedPositions( j ) );

                    rbfCoarse->interpolate2( valuesCoarse, valuesInterpolationCoarse );

                    double epsilon = std::sqrt( SMALL );
                    double error = ( valuesInterpolationCoarse.array() - this->values.array() ).matrix().norm() / (this->values.norm() + epsilon);
                    bool convergence = error < tolLivePointSelection;

                    if ( convergence )
                        reselection = false;

                    Info << "RBF interpolation coarsening: error = " << error << ", tol = " << tolLivePointSelection << ", reselection = ";

                    if ( reselection )
                        Info << "true";
                    else
                        Info << "false";

                    Info << endl;
                }

                if ( reselection )
                {
                    greedySelection( this->values );

                    rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
                }
            }
            else
            if ( !rbf->computed )
            {
                // Unit displacement of control points
                matrix unitDisplacement( positions.rows(), positions.cols() );
                unitDisplacement.setZero();

                assert( unitDisplacement.rows() >= nbMovingFaceCenters );

                if ( nbMovingFaceCenters == 0 )
                    unitDisplacement.fill( 1 );
                else
                    for ( int i = 0; i < nbMovingFaceCenters; i++ )
                        for ( int j = 0; j < unitDisplacement.cols(); j++ )
                            unitDisplacement( i, j ) = 1;

                greedySelection( unitDisplacement );

                rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
            }

            rbf::matrix selectedValues( selectedPositions.rows(), values.cols() );

            for ( int i = 0; i < selectedValues.rows(); i++ )
                selectedValues.row( i ) = values.row( selectedPositions( i ) );

            usedValues = selectedValues;
        }
        else
        {
            if ( !rbf->computed )
            {
                rbf->compute( positions, positionsInterpolation );
                rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
            }
        }

        usedValues.conservativeResize( usedValues.rows() - nbStaticFaceCentersRemove, usedValues.cols() );

        rbf->interpolate( usedValues, valuesInterpolation );
    }

    void RBFCoarsening::setNbMovingAndStaticFaceCenters(
        int nbMovingFaceCenters,
        int nbStaticFaceCenters
        )
    {
        nbStaticFaceCentersRemove = nbStaticFaceCenters;
        this->nbMovingFaceCenters = nbMovingFaceCenters;

        if ( enabled )
        {
            // Determine the number of selected static face centers
            nbStaticFaceCentersRemove = 0;

            for ( int i = 0; i < selectedPositions.rows(); i++ )
                if ( selectedPositions( i ) >= nbMovingFaceCenters )
                    nbStaticFaceCentersRemove++;
        }
    }
}
