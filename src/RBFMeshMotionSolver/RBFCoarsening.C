
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */
#include <sys/stat.h>
#include <iomanip>
#include "RBFCoarsening.H"

namespace rbf
{
    Foam::debug::debugSwitch RBFCoarsening::debug( "RBFInterpolation", 0 );

    RBFCoarsening::RBFCoarsening( std::shared_ptr<RBFInterpolation> rbf )
        :
        rbf( rbf ),
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction, rbf->polynomialTerm, rbf->cpu ) ) ),
        enabled( false ),
        livePointSelection( false ),
        livePointSelectionSumValues( false ),
        tol( 0 ),
        tolLivePointSelection( 0 ),
        coarseningMinPoints( 0 ),
        coarseningMaxPoints( 0 ),
        twoPointSelection( false ),
        surfaceCorrection( false ),
        surfaceCorrectionRbfFunction(),
        ratioRadiusError( 10 ),
        minCorrectionRadius( 1.0 ),
        surfaceCorrectionRadius( -1.0 ),
        cleanReselection( true ),
        exportTxt( false ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        errorInterpolationCoarse(),
        closestBoundaryIndexCorrection(),
        valuesCorrection(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 ),
        computed( false )
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
        bool exportTxt
        )
        :
        rbf( rbf ),
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction, rbf->polynomialTerm, rbf->cpu ) ) ),
        enabled( enabled ),
        livePointSelection( livePointSelection ),
        livePointSelectionSumValues( livePointSelectionSumValues ),
        tol( tol ),
        tolLivePointSelection( tolLivePointSelection ),
        coarseningMinPoints( coarseningMinPoints ),
        coarseningMaxPoints( coarseningMaxPoints ),
        twoPointSelection( false ),
        surfaceCorrection( false ),
        surfaceCorrectionRbfFunction(),
        ratioRadiusError( 10.0 ),
        surfaceCorrectionRadius( -1.0 ),
        cleanReselection( true ),
        exportTxt( exportTxt ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        errorInterpolationCoarse(),
        closestBoundaryIndexCorrection(),
        valuesCorrection(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 ),
        computed( false )
    {
        assert( rbf );
        assert( coarseningMinPoints <= coarseningMaxPoints );
        assert( coarseningMinPoints > 0 );
        assert( coarseningMaxPoints > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( tolLivePointSelection > 0 );
        assert( tolLivePointSelection < 1 );

        // If unit displacement do not use polynomial for selection
        if ( enabled && !livePointSelection && rbf->polynomialTerm )
        {
            WarningIn( "RBFCoarsening::RBFCoarsening" )
            << "Unit displacement is combined with polynomial addition into RBF interpolation. Could cause 'strange' results." << endl;
        }
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
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction, rbf->polynomialTerm, rbf->cpu ) ) ),
        enabled( enabled ),
        livePointSelection( livePointSelection ),
        livePointSelectionSumValues( livePointSelectionSumValues ),
        tol( tol ),
        tolLivePointSelection( tolLivePointSelection ),
        coarseningMinPoints( coarseningMinPoints ),
        coarseningMaxPoints( coarseningMaxPoints ),
        twoPointSelection( twoPointSelection ),
        surfaceCorrection( false ),
        surfaceCorrectionRbfFunction(),
        ratioRadiusError( 10.0 ),
        minCorrectionRadius( 1.0 ),
        surfaceCorrectionRadius( -1.0 ),
        cleanReselection( true ),
        exportTxt( exportTxt ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        errorInterpolationCoarse(),
        closestBoundaryIndexCorrection(),
        valuesCorrection(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 ),
        computed( false )
    {
        assert( rbf );
        assert( coarseningMinPoints <= coarseningMaxPoints );
        assert( coarseningMinPoints > 0 );
        assert( coarseningMaxPoints > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( tolLivePointSelection > 0 );
        assert( tolLivePointSelection < 1 );

        // If unit displacement do not use polynomial for selection
        if ( enabled && !livePointSelection && rbf->polynomialTerm )
        {
            WarningIn( "RBFCoarsening::RBFCoarsening" )
            << "Unit displacement is combined with polynomial addition into RBF interpolation. Could cause 'strange' results." << endl;
        }
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
        bool surfaceCorrection,
        std::shared_ptr<RBFFunctionInterface> surfaceCorrectionRbfFunction,
        double ratioRadiusError,
        double minCorrectionRadius,
        bool exportTxt
        )
        :
        rbf( rbf ),
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction, rbf->polynomialTerm, rbf->cpu ) ) ),
        enabled( enabled ),
        livePointSelection( livePointSelection ),
        livePointSelectionSumValues( livePointSelectionSumValues ),
        tol( tol ),
        tolLivePointSelection( tolLivePointSelection ),
        coarseningMinPoints( coarseningMinPoints ),
        coarseningMaxPoints( coarseningMaxPoints ),
        twoPointSelection( twoPointSelection ),
        surfaceCorrection( surfaceCorrection ),
        surfaceCorrectionRbfFunction( surfaceCorrectionRbfFunction ),
        ratioRadiusError( ratioRadiusError ),
        minCorrectionRadius( minCorrectionRadius ),
        surfaceCorrectionRadius( -1.0 ),
        cleanReselection( true ),
        exportTxt( exportTxt ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        errorInterpolationCoarse(),
        closestBoundaryIndexCorrection(),
        valuesCorrection(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 ),
        computed( false )
    {
        assert( rbf );
        assert( coarseningMinPoints <= coarseningMaxPoints );
        assert( coarseningMinPoints > 0 );
        assert( coarseningMaxPoints > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( tolLivePointSelection > 0 );
        assert( tolLivePointSelection < 1 );

        // If unit displacement do not use polynomial for selection
        if ( enabled && !livePointSelection && rbf->polynomialTerm )
        {
            WarningIn( "RBFCoarsening::RBFCoarsening" )
            << "Unit displacement is combined with polynomial addition into RBF interpolation. Could cause 'strange' results." << endl;
        }
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
        bool surfaceCorrection,
        std::shared_ptr<RBFFunctionInterface> surfaceCorrectionRbfFunction,
        double ratioRadiusError,
        double minCorrectionRadius,
        double surfaceCorrectionRadius,
        bool cleanReselection,
        bool exportTxt
        )
        :
        rbf( rbf ),
        rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction, rbf->polynomialTerm, rbf->cpu ) ) ),
        enabled( enabled ),
        livePointSelection( livePointSelection ),
        livePointSelectionSumValues( livePointSelectionSumValues ),
        tol( tol ),
        tolLivePointSelection( tolLivePointSelection ),
        coarseningMinPoints( coarseningMinPoints ),
        coarseningMaxPoints( coarseningMaxPoints ),
        twoPointSelection( twoPointSelection ),
        surfaceCorrection( surfaceCorrection ),
        surfaceCorrectionRbfFunction( surfaceCorrectionRbfFunction ),
        ratioRadiusError( ratioRadiusError ),
        minCorrectionRadius( minCorrectionRadius ),
        surfaceCorrectionRadius( surfaceCorrectionRadius ),
        cleanReselection( cleanReselection ),
        exportTxt( exportTxt ),
        selectedPositions(),
        nbStaticFaceCentersRemove( 0 ),
        positions(),
        positionsInterpolation(),
        values(),
        errorInterpolationCoarse(),
        closestBoundaryIndexCorrection(),
        valuesCorrection(),
        nbMovingFaceCenters( 0 ),
        fileExportIndex( 0 ),
        computed( false )
    {
        assert( rbf );
        assert( coarseningMinPoints <= coarseningMaxPoints );
        assert( coarseningMinPoints > 0 );
        assert( coarseningMaxPoints > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( tolLivePointSelection > 0 );
        assert( tolLivePointSelection < 1 );

        // If unit displacement do not use polynomial for selection
        if ( enabled && !livePointSelection && rbf->polynomialTerm )
        {
            WarningIn( "RBFCoarsening::RBFCoarsening" )
            << "Unit displacement is combined with polynomial addition into RBF interpolation. Could cause 'strange' results." << endl;
        }

        if ( surfaceCorrection && surfaceCorrectionRadius > 0 )
        {
            WarningIn( "RBFCoarsening::RBFCoarsening" )
            << "Using surfaceCorrection and radius is manually set by surfaceCorrectionRadius: " << surfaceCorrectionRadius << ". This could lead to bad/invalid meshes." << endl;
        }

        if ( debug == 3 )
        {
            std::string filename = "totalInterpolationTimings.txt";
            std::ofstream timingFile( filename, std::ofstream::app );

            if ( Pstream::myProcNo() == 0 && not rbf->polynomialTerm )
                timingFile << std::setprecision( 9 ) << "t_total, t_error, t_reselect, t_interpolate, t_correct" << "\n";
        }
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
        assert( !computed );

        this->positions = positions;
        this->positionsInterpolation = positionsInterpolation;

        computed = true;
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
            int maxNbPoints = std::min( coarseningMaxPoints, static_cast<int>( positions.rows() ) );
            int minPoints = std::min( coarseningMinPoints, static_cast<int>( positions.rows() ) );
            double error = 0;
            double errorMax = 0;
            rbf::vector errorList( positions.rows() );

            // Clean selected points if requested
            if ( cleanReselection || selectedPositions.rows() == maxNbPoints || selectedPositions.rows() == 0 )
            {
                selectedPositions.resize( 2 );
                if( livePointSelection )
                {
                    // Select the point with the largest displacment
                    int maxDisplacementIndex = -1;
                    ( values.rowwise().norm() ).maxCoeff( &maxDisplacementIndex );
                    //Add first point
                    selectedPositions( 0 ) = maxDisplacementIndex;
                }
                else
                {
                    //With unit displacement, first point is point with largest radius from origin and displacement > 0
                    vector rad = positions.rowwise().norm();

                    int maxRadiusFromOriginIndex = -1;
                    double maxRadius = -1;
                    for(int i=0;i<rad.rows();i++)
                    {
                        if(rad(i) > maxRadius &&  values.row(i).norm() > SMALL )
                        {
                            maxRadius = rad(i);
                            maxRadiusFromOriginIndex = i;
                        }
                    }
                    selectedPositions( 0 ) = maxRadiusFromOriginIndex;
                }

                //Find point with largest distance from first point
                vector rad = ( positions - ( matrix::Constant( positions.rows(), 1, 1.0 )*positions.row( selectedPositions(0) ) ) ).rowwise().norm();
                int maxRadiusIndex = -1;
                double maxRadius = -1;
                for(int i=0;i<rad.rows();i++)
                {
                    if(rad(i) > maxRadius && (rad(i) < 1.0 - SMALL || rad(i) > 1.0 + SMALL))
                    {
                        maxRadius = rad(i);
                        maxRadiusIndex = i;
                    }
                }
                //Add second point
                selectedPositions( 1 ) = maxRadiusIndex;

            }

            assert( positions.rows() >= selectedPositions.rows() );

            rbf::matrix positionsInterpolationCoarse = positions;

            // Create RBF interpolator

            // Run the greedy algorithm
            double runTimeInterpolate = 0.0;
            double runTimeError = 0.0;
            double runTimeConvergence = 0.0;
            bool addedSecondPoint = false;
            int counter = selectedPositions.rows();
            int counter0 = counter;

            while ( true )
            {
                // Build the matrices used for the RBF interpolation
                rbf::matrix positionsCoarse( counter, positions.cols() );
                rbf::matrix valuesCoarse( positionsCoarse.rows(), positionsCoarse.cols() );
                rbf::matrix valuesInterpolationCoarse( positionsInterpolationCoarse.rows(), positionsInterpolationCoarse.cols() );

                for ( int j = 0; j < selectedPositions.rows(); j++ )
                {
                    positionsCoarse.row( j ) = positions.row( selectedPositions( j ) );
                    valuesCoarse.row( j ) = values.row( selectedPositions( j ) );
                }

                std::clock_t t = std::clock();

                // Perform the RBF interpolation.
                rbfCoarse->interpolate( positionsCoarse, positionsInterpolationCoarse, valuesCoarse, valuesInterpolationCoarse );

                if ( debug > 0 )
                {
                    t = std::clock() - t;
                    runTimeInterpolate += static_cast<float>(t) / CLOCKS_PER_SEC;
                    t = std::clock();
                }

                // Evaluate the error
                for ( int j = 0; j < valuesInterpolationCoarse.rows(); j++ )
                    errorList( j ) = ( valuesInterpolationCoarse.row( j ) - values.row( j ) ).norm();

                // Select the point with the largest error which is not already selected.
                int index = -1;
                double largestError = errorList.maxCoeff( &index );

                // Additional function to check whether the largestError = 0 (<SMALL) and do select next consecutive point
                if ( largestError < SMALL )
                {
                    index = selectedPositions.rows();
                }

                int index2 = -1;
                double largestError2 = -1;

                // selected point with largest error in opposite direction (more than 90 degrees differenc in direction)
                if ( twoPointSelection )
                {
                    vector largestErrorVector = valuesInterpolationCoarse.row( index ) - values.row( index );

                    for ( int j = 0; j < errorList.rows(); j++ )
                    {
                        vector errorVector = valuesInterpolationCoarse.row( j ) - values.row( j );

                        if ( largestErrorVector.dot( errorVector ) < -SMALL && largestError2 < errorList( j ) )
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
                }

                double epsilon = std::sqrt( SMALL );
                error = (errorList).norm() / (values.norm() + epsilon);
                errorMax = largestError / ( ( values.rowwise().norm() ).maxCoeff() + epsilon );

                // bool convergence = (error < tol && counter >= minPoints) || counter >= maxNbPoints;//only take 2-norm
                bool convergence = (error < tol && errorMax < tol && counter >= minPoints) || counter >= maxNbPoints;//take both 2-norm and Inf-norm

                if ( convergence )
                {
                    if ( livePointSelection )
                        errorInterpolationCoarse = valuesInterpolationCoarse - values;

                    break;
                }

                // Break if maximum point are reached
                if ( counter >= maxNbPoints )
                {
                    if ( livePointSelection )
                        errorInterpolationCoarse = valuesInterpolationCoarse - values;

                    break;
                }

                selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
                selectedPositions( selectedPositions.rows() - 1 ) = index;
                counter++;

                // Add second point if possible
                if ( twoPointSelection && index2 >= 0 && index != index2 )
                {
                    addedSecondPoint = true;
                    selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
                    selectedPositions( selectedPositions.rows() - 1 ) = index2;
                    counter++;
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
                Info << "RBFCoarsening::debug 2. find largest error = " << runTimeError << " s" << ". Added second point = " << addedSecondPoint << endl;
                Info << "RBFCoarsening::debug 3. convergence check = " << runTimeConvergence << " s" << endl;

                if ( livePointSelection )
                {
                    Info << "RBFCoarsening::debug 4. absolute max error = " << ( errorInterpolationCoarse.rowwise().norm() ).maxCoeff() << " m" << endl;
                }

                if ( !cleanReselection )
                    Info << "RBFCoarsening::debug 5. re-used selection and added " << counter - counter0 << " points to control points." << endl;
            }

            Info << "RBF interpolation coarsening: selected " << selectedPositions.rows() << "/" << positions.rows() << " points, 2-norm(error) = "
                 << error << ", max(error) = " << errorMax << ", tol = " << tol << endl;

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
        std::clock_t t = std::clock();
        std::clock_t tp = std::clock();
        double runTimeError = 0.0;
        double runTimeReselect = 0.0;
        double runTimeInterpolate = 0.0;
        double runTimeCorrect = 0.0;

        matrix usedValues = values;

        if ( enabled )
        {
            if ( livePointSelection )
            {
                if ( debug == 3 )
                {
                    tp = std::clock();
                }

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

                    errorInterpolationCoarse = valuesInterpolationCoarse - this->values;
                    double error = (errorInterpolationCoarse).matrix().norm() / (this->values.norm() + epsilon);
                    double errorMax = ( errorInterpolationCoarse.rowwise().norm() ).maxCoeff() / ( ( this->values.rowwise().norm() ).maxCoeff() + epsilon );

                    // bool convergence = error < tolLivePointSelection;//only take the 2-norm error
                    bool convergence = (error < tolLivePointSelection && errorMax < tolLivePointSelection);//take both 2-norm and Inf-norm error as convergence

                    if ( convergence )
                        reselection = false;

                    Info << "RBF interpolation coarsening: 2-norm(error) = " << error << ", max(error) = " << errorMax << ", tol = " << tolLivePointSelection << ", reselection = ";

                    if ( reselection )
                        Info << "true";
                    else
                        Info << "false";

                    Info << endl;

                    // If debug is 2: Print out surface error to file
                    if ( debug >= 2 && Pstream::myProcNo() == 0 && not rbf->polynomialTerm )
                    {
                        Info << "Debug 2: 2-norm(error|disp) = " << (errorInterpolationCoarse).matrix().norm() << "|" << this->values.norm() << ", Inf-norm(error|disp) = " << (errorInterpolationCoarse.rowwise().norm() ).maxCoeff() << "|" << (this->values.rowwise().norm() ).maxCoeff() << endl;

                        std::string filename = "liveSelection-rbf-surfaceError.txt";

                        std::ofstream surfaceErrorFile( filename, std::ofstream::app );

                        if ( surfaceErrorFile.is_open() )
                        {
                            surfaceErrorFile << error << ", " << errorMax;

                            if ( !reselection )
                            {
                                surfaceErrorFile << ", " << valuesCoarse.rows() << ", " << valuesInterpolationCoarse.rows() << ", " << error << ", " << errorMax << ", " << reselection << ", " << (errorInterpolationCoarse).matrix().norm() << ", " << this->values.norm() << ", " << (errorInterpolationCoarse.rowwise().norm() ).maxCoeff() << ", " << (this->values.rowwise().norm() ).maxCoeff() << "\n";
                            }
                        }

                        surfaceErrorFile.close();
                    }
                }

                if ( debug == 3 )
                {
                    tp = std::clock() - tp;
                    runTimeError = static_cast<float>(tp) / CLOCKS_PER_SEC;
                    tp = std::clock();
                }

                if ( reselection )
                {
                    greedySelection( this->values );

                    rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );

                    // If debug is 2: Print out surface error to file
                    if ( debug >= 2 && Pstream::myProcNo() == 0 && not rbf->polynomialTerm )
                    {
                        std::string filename = "liveSelection-rbf-surfaceError.txt";
                        struct stat buffer;
                        bool fileExists = (stat( filename.c_str(), &buffer ) == 0);

                        // ==== START Re-calculate error ==== //
                        rbf::matrix valuesCoarse( selectedPositions.rows(), this->values.cols() );
                        rbf::matrix valuesInterpolationCoarse( positions.rows(), valuesInterpolation.cols() );

                        for ( int j = 0; j < selectedPositions.rows(); j++ )
                            valuesCoarse.row( j ) = this->values.row( selectedPositions( j ) );

                        rbfCoarse->interpolate2( valuesCoarse, valuesInterpolationCoarse );

                        double epsilon = std::sqrt( SMALL );

                        errorInterpolationCoarse = valuesInterpolationCoarse - this->values;
                        double error = (errorInterpolationCoarse).matrix().norm() / (this->values.norm() + epsilon);
                        double errorMax = ( errorInterpolationCoarse.rowwise().norm() ).maxCoeff() / ( ( this->values.rowwise().norm() ).maxCoeff() + epsilon );

                        // ==== END Re-calculate error ==== //

                        if ( !fileExists )
                        {
                            std::ofstream surfaceErrorFile( filename, std::ofstream::app );

                            surfaceErrorFile << "eps_2, eps_max, Nc, Nb, eps_2, eps_max, reselected, error_2, dx_2, error_max, dx_max" << "\n";
                            if ( surfaceErrorFile.is_open() )
                            {
                                surfaceErrorFile << error << ", " << errorMax << ", " << valuesCoarse.rows() << ", " << valuesInterpolationCoarse.rows() << ", " << error << ", " << errorMax << ", " << reselection << ", " << (errorInterpolationCoarse).matrix().norm() << ", " << this->values.norm() << ", " << (errorInterpolationCoarse.rowwise().norm() ).maxCoeff() << ", " << (this->values.rowwise().norm() ).maxCoeff() << "\n";
                            }

                            surfaceErrorFile.close();
                        }
                        else
                        {
                            std::ofstream surfaceErrorFile( filename, std::ofstream::app );

                            if ( surfaceErrorFile.is_open() )
                            {
                                surfaceErrorFile << ", " << valuesCoarse.rows() << ", " << valuesInterpolationCoarse.rows() << ", " << error << ", " << errorMax << ", " << reselection << ", " << (errorInterpolationCoarse).matrix().norm() << ", " << this->values.norm() << ", " << (errorInterpolationCoarse.rowwise().norm() ).maxCoeff() << ", " << (this->values.rowwise().norm() ).maxCoeff() << "\n";
                            }
                        }
                    }
                }

                if ( debug == 3 )
                {
                    tp = std::clock() - tp;
                    runTimeReselect = static_cast<float>(tp) / CLOCKS_PER_SEC;
                    tp = std::clock();
                }
            }
            else
            if ( !rbf->computed )//unit displacement
            {
                if ( debug == 3 )
                {
                    tp = std::clock();
                }

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

                if ( debug == 3 )
                {
                    tp = std::clock() - tp;
                    runTimeReselect = static_cast<float>(tp) / CLOCKS_PER_SEC;
                    tp = std::clock();
                }

                if( surfaceCorrection )
                {
                    if ( livePointSelectionSumValues )
                    {
                        if ( this->values.cols() != values.cols() )
                            this->values = values;
                        else
                            this->values.array() += values.array();
                    }

                    // Construct values to interpolate based on unit displacement selected points
                    rbf::matrix valuesCoarse( selectedPositions.rows(), this->values.cols() );
                    rbf::matrix valuesInterpolationCoarse( positions.rows(), valuesInterpolation.cols() );

                    for ( int j = 0; j < selectedPositions.rows(); j++ )
                        valuesCoarse.row( j ) = this->values.row( selectedPositions( j ) );

                    // This will return the displaced surface in valuesInterpolationCoarse
                    rbfCoarse->interpolate2( valuesCoarse, valuesInterpolationCoarse );

                    //Calculate error for surfaceCorrection
                    errorInterpolationCoarse = valuesInterpolationCoarse - this->values;

                    if ( debug > 0 )
                    {
                        rbf::vector errorList( errorInterpolationCoarse.rows() );

                        // Evaluate the error
                        for ( int j = 0; j < errorInterpolationCoarse.rows(); j++ )
                            errorList( j ) = ( errorInterpolationCoarse.row( j ) ).norm();

                        double epsilon = std::sqrt( SMALL );
                        double error = (errorList).norm() / (this->values.norm() + epsilon);
                        double errorMax = errorList.maxCoeff() / ( ( this->values.rowwise().norm() ).maxCoeff() + epsilon );

                        Info << "RBFCoarsening::UnitDisplacement::debug 1: " << "2-norm error = " << error << ", max error = " << errorMax << endl;

                        // If debug is 2: Print out surface error to file
                        if ( debug >= 2 && Pstream::myProcNo() == 0 && not rbf->polynomialTerm )
                        {
                            std::string filename = "unitSelection-rbf-surfaceError.txt";
                            std::ofstream surfaceErrorFile( filename, std::ofstream::app );

                            if ( surfaceErrorFile.is_open() )
                            {
                                surfaceErrorFile << error << ", " << errorMax << ", " << valuesCoarse.rows() << ", " << valuesInterpolationCoarse.rows() << "\n";
                            }

                            surfaceErrorFile.close();
                        }
                    }
                }

                rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );

                if ( debug == 3 )
                {
                    tp = std::clock() - tp;
                    runTimeError = static_cast<float>(tp) / CLOCKS_PER_SEC;
                    tp = std::clock();
                }
            }
            else // This means there is unit displacement used, but rbf is already computed. Only used for debug things to track error.
            {
                if ( debug == 3 )
                {
                    tp = std::clock();
                }

                if ( surfaceCorrection )
                {
                    if ( livePointSelectionSumValues )
                    {
                        if ( this->values.cols() != values.cols() )
                            this->values = values;
                        else
                            this->values.array() += values.array();
                    }

                    // Construct values to interpolate based on unit displacement selected points
                    rbf::matrix valuesCoarse( selectedPositions.rows(), values.cols() );
                    rbf::matrix valuesInterpolationCoarse( positions.rows(), valuesInterpolation.cols() );

                    for ( int j = 0; j < selectedPositions.rows(); j++ )
                        valuesCoarse.row( j ) = this->values.row( selectedPositions( j ) );

                    // This will return the displaced surface in valuesInterpolationCoarse
                    rbfCoarse->interpolate2( valuesCoarse, valuesInterpolationCoarse );

                    //Calculate error for surfaceCorrection
                    errorInterpolationCoarse = valuesInterpolationCoarse - this->values;

                    if ( debug > 0 )
                    {
                        // Evaluate the error
                        rbf::vector errorList( positions.rows() );
                        for ( int j = 0; j < errorInterpolationCoarse.rows(); j++ )
                            errorList( j ) = ( errorInterpolationCoarse.row( j ) ).norm();

                        double epsilon = std::sqrt( SMALL );
                        double error = (errorList).norm() / (this->values.norm() + epsilon);
                        double errorMax = errorList.maxCoeff() / ( ( this->values.rowwise().norm() ).maxCoeff() + epsilon );

                        Info << "RBFCoarsening::UnitDisplacement::debug 1: " << "2-norm error = " << error << ", max error = " << errorMax << endl;

                        // If debug is 2: Print out surface error to file
                        if ( debug >= 2 && Pstream::myProcNo() == 0 && not rbf->polynomialTerm )
                        {
                            std::string filename = "unitSelection-rbf-surfaceError.txt";
                            std::ofstream surfaceErrorFile( filename, std::ofstream::app );

                            if ( surfaceErrorFile.is_open() )
                            {
                                surfaceErrorFile << error << ", " << errorMax << ", " << valuesCoarse.rows() << ", " << valuesInterpolationCoarse.rows() << "\n";
                            }

                            surfaceErrorFile.close();
                        }
                    }
                }

                if ( debug == 3 )
                {
                    tp = std::clock() - tp;
                    runTimeReselect = static_cast<float>(tp) / CLOCKS_PER_SEC;
                    tp = std::clock();
                }
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

        if ( debug == 3 )
        {
            tp = std::clock() - tp;
            runTimeInterpolate = static_cast<float>(tp) / CLOCKS_PER_SEC;
            tp = std::clock();
        }

        // start doing correction of surface is requested
        if ( surfaceCorrection )
        {
            correctSurface( valuesInterpolation, errorInterpolationCoarse );
        }

        //DEBUG STATEMENT
        if ( debug == 3 )
        {
            tp = std::clock() - tp;
            runTimeCorrect = static_cast<float>(tp) / CLOCKS_PER_SEC;
            tp = std::clock();
        }
        //DEBUG STATEMENT
        if ( debug == 3 )
        {
            t = std::clock() - t;
            double runTimeINTP = static_cast<float>(t) / CLOCKS_PER_SEC;
            t = std::clock();
            Info << "RBFCoarsening::interpolate::debug 1. total time = " << runTimeINTP << " s" << endl;

            // write to file
            if ( Pstream::myProcNo() == 0 && not rbf->polynomialTerm )
            {
                std::string filename = "totalInterpolationTimings.txt";
                std::ofstream timingFile( filename, std::ofstream::app );

                if ( timingFile.is_open() )
                {
                    timingFile << std::setprecision( 9 ) << runTimeINTP << ", " << runTimeError << ", " << runTimeReselect << ", " << runTimeInterpolate << ", " << runTimeCorrect << "\n";
                }

                timingFile.close();
            }
        }
    }

    void RBFCoarsening::correctSurface( matrix & valuesInterpolation, const matrix & surfaceError )
    {
        double maxError = ( surfaceError.rowwise().norm() ).maxCoeff();
        //ensure that only performed when there is an error bigger than 0
        if( maxError > SMALL )
        {
            if ( valuesCorrection.rows() == 0 )
            {
                valuesCorrection.conservativeResize( valuesInterpolation.rows(), valuesInterpolation.cols() );
                valuesCorrection.setZero();
            }

            double R = 1.0;

            if ( surfaceCorrectionRadius > 0 )
            {
                R = surfaceCorrectionRadius;
            }
            else
            {
                double Rerror = ratioRadiusError * maxError;
                double Rwall = ratioRadiusError * minCorrectionRadius;
                R = max( Rerror, Rwall );
            }

            if ( debug > 0 )
            {
                if ( surfaceCorrectionRadius > 0 )
                {
                    Info << nl << "RBFCoarsening::correctSurface::debug 0: R = " << R << endl;
                }
                else
                {
                    Info << nl << "RBFCoarsening::correctSurface::debug 0: ratioRadiusError = " << ratioRadiusError << ", R = " << R << endl;
                }
            }

            // Find nearest boundary point for each internal point. Do this only the first time
            vector closestBoundaryRadius( positionsInterpolation.rows() );

            if ( closestBoundaryIndexCorrection.rows() == 0 )
            {
                closestBoundaryIndexCorrection.conservativeResize( positionsInterpolation.rows() );
                std::clock_t t = std::clock();
                double runTimeNN = 0;

                for ( int i = 0; i < positionsInterpolation.rows(); i++ )
                {
                    double smallestRadius = GREAT;
                    int boundaryIndex = -1;

                    for ( int j = 0; j < positions.rows(); j++ )
                    {
                        double radius = ( positions.row( j ) - positionsInterpolation.row( i ) ).norm();

                        if ( radius < smallestRadius )
                        {
                            boundaryIndex = j;
                            smallestRadius = radius;
                        }
                    }

                    closestBoundaryIndexCorrection( i ) = boundaryIndex;
                    closestBoundaryRadius( i ) = smallestRadius;
                }

                if ( debug > 0 )
                {
                    t = std::clock() - t;
                    runTimeNN += static_cast<float>(t) / CLOCKS_PER_SEC;
                    t = std::clock();
                    Info << "RBFCoarsening::correctSurface::debug 1. nearest neighbour selection = " << runTimeNN << " s" << endl;
                }
            }
            else
            {
                for ( int i = 0; i < positionsInterpolation.rows(); i++ )
                {
                    closestBoundaryRadius( i ) = ( positions.row( closestBoundaryIndexCorrection( i ) ) - positionsInterpolation.row( i ) ).norm();
                }
            }

            // Start doing the correction
            std::clock_t t = std::clock();
            double runTimeCorr = 0;
            surfaceCorrectionRbfFunction->setRadius( R );

            for ( int i = 0; i < positionsInterpolation.rows(); i++ )
            {
                matrix fEval = -( surfaceCorrectionRbfFunction->evaluate( closestBoundaryRadius( i ) ) ) * surfaceError.row( closestBoundaryIndexCorrection( i ) );
                valuesInterpolation.row( i ) += ( fEval - valuesCorrection.row( i ) );
                valuesCorrection.row( i ) = fEval;
            }

            if ( debug > 0 )
            {
                t = std::clock() - t;
                runTimeCorr += static_cast<float>(t) / CLOCKS_PER_SEC;
                t = std::clock();
                Info << "RBFCoarsening::correctSurface::debug 2. correction evaluation = " << runTimeCorr << " s" << endl;
            }
        }
        else
        {
            Info << "RBF Coarsening: No surface correction is performed since error on boundary is zero." << nl << endl;
        }
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
