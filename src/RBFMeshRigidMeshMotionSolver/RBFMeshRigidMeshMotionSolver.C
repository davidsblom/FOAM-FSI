
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFMeshRigidMeshMotionSolver.H"

using namespace Foam;

defineTypeNameAndDebug( RBFMeshRigidMeshMotionSolver, 0 );

addToRunTimeSelectionTable
(
    motionSolver,
    RBFMeshRigidMeshMotionSolver,
    dictionary
);

RBFMeshRigidMeshMotionSolver::RBFMeshRigidMeshMotionSolver(
    const polyMesh & mesh,
    Istream & msData
    )
    :
    RBFMeshMotionSolver( mesh, msData ),
    translationAmplitude( 0 ),
    translationFrequency( 0 ),
    translationDirection( Foam::vector::zero ),
    rotationAmplitude( 0 ),
    rotationFrequency( 0 ),
    rotationOrigin( Foam::vector::zero ),
    oldTransformation( Foam::vector::zero ),
    told( 0 )
{
    dictionary & rigidBodyMotionDict = subDict( "rigidBodyMotion" );
    translationAmplitude = readScalar( rigidBodyMotionDict.lookup( "translationAmplitude" ) );
    translationFrequency = readScalar( rigidBodyMotionDict.lookup( "translationFrequency" ) );
    rotationAmplitude = readScalar( rigidBodyMotionDict.lookup( "rotationAmplitude" ) );
    rotationFrequency = readScalar( rigidBodyMotionDict.lookup( "rotationFrequency" ) );
    translationDirection = rigidBodyMotionDict.lookup( "translationDirection" );
    rotationOrigin = rigidBodyMotionDict.lookup( "rotationOrigin" );

    string filename = static_cast<std::string>( mesh.time().rootPath() ) + "/" + static_cast<std::string>( mesh.time().caseName() ) + "/constant/fsi.yaml";

    YAML::Node config = YAML::LoadFile( filename );

    assert( config["sdc"] || config["esdirk"] );

    if ( config["sdc"] )
    {
        YAML::Node sdcConfig( config["sdc"] );
        assert( sdcConfig["convergence-tolerance"] );
        assert( sdcConfig["number-of-points"] );
        assert( sdcConfig["quadrature-rule"] );

        int n = sdcConfig["number-of-points"].as<int>();
        scalar tol = sdcConfig["convergence-tolerance"].as<scalar>();
        std::string quadratureRule = sdcConfig["quadrature-rule"].as<std::string>();

        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( quadratureRule, n, tol ) );
    }

    if ( config["esdirk"] )
    {
        YAML::Node esdirkConfig( config["esdirk"] );

        assert( esdirkConfig["method"] );

        std::string method = esdirkConfig["method"].as<std::string>();

        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::ESDIRK( method ) );
    }

    int nbTimeSteps = mesh.time().endTime().value() / mesh.time().deltaT().value();
    scalar dt = mesh.time().deltaT().value();
    fsi::vector q0 ( 2 );
    q0 << 1, 0;
    scalar amplitude = 0;
    scalar frequency = 1;
    scalar m = 1;
    scalar k = 100;
    oscillator = std::shared_ptr<Oscillator>( new Oscillator( nbTimeSteps, dt, q0, amplitude, frequency, m, k ) );
}

RBFMeshRigidMeshMotionSolver::~RBFMeshRigidMeshMotionSolver()
{}

Foam::vector RBFMeshRigidMeshMotionSolver::calcVelocity()
{
    assert( timeIntegrationScheme );

    scalar t = mesh().time().value();
    scalar dt = mesh().time().deltaT().value();

    fsi::vector rhs( 2 ), result( 2 ), f( 2 ), qold;

    if ( not corrector && k == 0 )
        timeIntegrationScheme->setOldSolution( mesh().time().timeIndex(), oscillator->sol );

    timeIntegrationScheme->getSourceTerm( corrector, k, dt, rhs, qold );

    oscillator->implicitSolve( corrector, k, k, t, dt, qold, rhs, f, result );

    timeIntegrationScheme->setFunction( k, f, result );

    scalar disp = oscillator->sol(0);
    // change reference frame
    disp *= -1;
    disp += 1;
    disp *= translationAmplitude;

    Foam::vector transformation = Foam::vector::zero;
    transformation[1] = disp - oldTransformation[1];

    Info << "dudt(1) = " << oscillator->sol(0) << endl;

    told = t;
    oldTransformation[1] = disp;

    return transformation;
}

void RBFMeshRigidMeshMotionSolver::solve()
{
    /*
     * RBF interpolator from face centers to local complete mesh vertices
     * The interpolation consists of the following steps:
     * 1. Build a matrix with the face center positions of the static patches and the moving patches
     * 2. Build a matrix with the positions of every vertex in the local mesh
     * 3. Build a matrix with the displacement/motion of the face center positions of the static patches and the moving patches
     * 4. Perform the interpolation from the face centers to the complete mesh
     * 5. Correct the mesh vertices of the static patches. Set these displacement to zero.
     * 6. Set the motion of the mesh vertices
     */

    /*
     * Step 1: Build a matrix with the face center positions of the static patches and the moving patches
     * The order of the matrix is defined as first a list of the moving patch face centers,
     * thereafter the static patch face centers. These are the control points used by the
     * radial basis function interpolation.
     * The control points should be exactly the same at each processor, and are therefore communicated
     * to each process. As only an absolute RBF interpolation is implemented, this communication step is only
     * performed once per simulation.
     * The global ordering of the data is first the information of the moving patches,
     * thereafter the static patches.
     */

    int nbFaceCenters = 0;
    int nbMovingFaceCenters = 0;
    int nbStaticFaceCenters = 0;
    int nbFixedFaceCenters = 0;

    std::vector<int> staticControlPointLabels;
    std::vector<int> staticControlGlobalPointLabels;
    std::vector<int> fixedControlPointLabels;
    std::vector<int> fixedControlGlobalPointLabels;
    std::vector<int> movingControlPointLabels;
    std::vector<int> movingControlPointPatchIds;
    std::vector<int> movingControlPointIndices;
    std::vector<int> movingControlGlobalPointLabels;

    labelList globalStaticPointsListEnabled( nbStaticFaceCenters, 0 );
    labelList globalFixedPointsListEnabled( nbFixedFaceCenters, 0 );
    labelList globalMovingPointsListEnabled( nbMovingFaceCenters, 0 );
    int globalStaticOffsetNonUnique = 0;
    int globalFixedOffsetNonUnique = 0;
    int globalMovingOffsetNonUnique = 0;

    if ( sum( nbGlobalFaceCenters ) == 0 )
    {
        // Determine the number of face centers
        // The total number of face centers is simply the sum of the face centers
        // on each processor.
        nbFaceCenters = 0;

        // First add the static patches, thereafter the fixed patches, and
        // the moving patches as last.

        forAll( staticPatchIDs, i )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[staticPatchIDs[i]].meshPoints();

            forAll( meshPoints, j )
            {
                if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                    continue;

                // Only add the static vertex point if it's not already added to the list
                std::vector<int>::iterator it;

                it = std::find( staticControlPointLabels.begin(), staticControlPointLabels.end(), meshPoints[j] );

                if ( it == staticControlPointLabels.end() || staticControlPointLabels.size() == 0 )
                {
                    staticControlPointLabels.push_back( meshPoints[j] );
                }
            }
        }

        nbStaticFaceCenters = staticControlPointLabels.size();

        forAll( fixedPatchIDs, i )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[fixedPatchIDs[i]].meshPoints();

            forAll( meshPoints, j )
            {
                if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                    continue;

                // Only add the static vertex point if it's not already added to the list
                // and if it's not a static point
                std::vector<int>::iterator it;

                it = std::find( staticControlPointLabels.begin(), staticControlPointLabels.end(), meshPoints[j] );

                if ( it == staticControlPointLabels.end() || staticControlPointLabels.size() == 0 )
                {
                    it = std::find( fixedControlPointLabels.begin(), fixedControlPointLabels.end(), meshPoints[j] );

                    if ( it == fixedControlPointLabels.end() || fixedControlPointLabels.size() == 0 )
                    {
                        fixedControlPointLabels.push_back( meshPoints[j] );
                    }
                }
            }
        }

        nbFixedFaceCenters = fixedControlPointLabels.size();

        if ( faceCellCenters )
        {
            forAll( movingPatchIDs, i )
            {
                nbMovingFaceCenters += mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres().size();
            }
        }

        if ( not faceCellCenters )
        {
            forAll( movingPatchIDs, patchI )
            {
                const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchI]].meshPoints();
                globalMovingPointsLabelList[movingPatchIDs[patchI]] = labelList( meshPoints.size(), 0 );

                forAll( meshPoints, j )
                {
                    if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                        continue;

                    // Only add the static vertex point if it's not already added to the list
                    // and if it's not a static point
                    std::vector<int>::iterator it;

                    it = std::find( staticControlPointLabels.begin(), staticControlPointLabels.end(), meshPoints[j] );

                    if ( it == staticControlPointLabels.end() || staticControlPointLabels.size() == 0 )
                    {
                        it = std::find( fixedControlPointLabels.begin(), fixedControlPointLabels.end(), meshPoints[j] );

                        if ( it == fixedControlPointLabels.end() || fixedControlPointLabels.size() == 0 )
                        {
                            it = std::find( movingControlPointLabels.begin(), movingControlPointLabels.end(), meshPoints[j] );

                            if ( it == movingControlPointLabels.end() || movingControlPointLabels.size() == 0 )
                            {
                                movingControlPointLabels.push_back( meshPoints[j] );
                                movingControlPointPatchIds.push_back( movingPatchIDs[patchI] );
                                movingControlPointIndices.push_back( j );
                                globalMovingPointsLabelList[movingPatchIDs[patchI]][j] = 1;
                            }
                        }
                    }
                }
            }

            nbMovingFaceCenters = movingControlPointLabels.size();
        }

        if ( Pstream::nProcs() == 1 )
        {
            globalStaticPointsListEnabled.resize( nbStaticFaceCenters );
            globalStaticPointsListEnabled = 1;
            globalFixedPointsListEnabled.resize( nbFixedFaceCenters );
            globalFixedPointsListEnabled = 1;
            globalMovingPointsListEnabled.resize( nbMovingFaceCenters );
            globalMovingPointsListEnabled = 1;
        }

        if ( Pstream::nProcs() > 1 )
        {
            assert( false );
            IOobject addrHeader
            (
                "pointProcAddressing",
                mesh().facesInstance(),
                mesh().meshSubDir,
                mesh(),
                IOobject::MUST_READ
            );

            assert( addrHeader.headerOk() );
            labelIOList pointProcAddressing( addrHeader );

            assert( pointProcAddressing.size() == mesh().points().size() );

            // Count the number of global static points including scalar points
            nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
            nbGlobalFixedFaceCenters[Pstream::myProcNo()] = nbFixedFaceCenters;
            nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
            reduce( nbGlobalStaticFaceCenters, sumOp<labelList>() );
            reduce( nbGlobalFixedFaceCenters, sumOp<labelList>() );
            reduce( nbGlobalMovingFaceCenters, sumOp<labelList>() );
            nbStaticFaceCenters = sum( nbGlobalStaticFaceCenters );
            nbFixedFaceCenters = sum( nbGlobalFixedFaceCenters );

            if ( not faceCellCenters )
                nbMovingFaceCenters = sum( nbGlobalMovingFaceCenters );

            // Construct a list with all the global point labels, thus including
            // also scalar points. Thereafter, construct a list of static control
            // list which indicates whether the point is already included or not.
            // Use this later to build a list of the unique static control points.
            labelList globalStaticPointsList( nbStaticFaceCenters, 0 );
            labelList globalFixedPointsList( nbFixedFaceCenters, 0 );
            labelList globalMovingPointsList( nbMovingFaceCenters, 0 );
            labelList globalMovingPointsPatchIds( nbMovingFaceCenters, 0 );
            labelList globalMovingPointsIndices( nbMovingFaceCenters, 0 );

            globalStaticOffsetNonUnique = 0;

            for ( int i = 0; i < Pstream::myProcNo(); i++ )
                globalStaticOffsetNonUnique += nbGlobalStaticFaceCenters[i];

            globalFixedOffsetNonUnique = 0;

            for ( int i = 0; i < Pstream::myProcNo(); i++ )
                globalFixedOffsetNonUnique += nbGlobalFixedFaceCenters[i];

            globalMovingOffsetNonUnique = 0;

            for ( int i = 0; i < Pstream::myProcNo(); i++ )
                globalMovingOffsetNonUnique += nbGlobalMovingFaceCenters[i];

            for ( unsigned int i = 0; i < staticControlPointLabels.size(); i++ )
            {
                globalStaticPointsList[i + globalStaticOffsetNonUnique] = pointProcAddressing[staticControlPointLabels[i]];
            }

            for ( unsigned int i = 0; i < fixedControlPointLabels.size(); i++ )
            {
                globalFixedPointsList[i + globalFixedOffsetNonUnique] = pointProcAddressing[fixedControlPointLabels[i]];
            }

            for ( unsigned int i = 0; i < movingControlPointLabels.size(); i++ )
            {
                globalMovingPointsList[i + globalMovingOffsetNonUnique] = pointProcAddressing[movingControlPointLabels[i]];
                globalMovingPointsPatchIds[i + globalMovingOffsetNonUnique] = movingControlPointPatchIds[i];
                globalMovingPointsIndices[i + globalMovingOffsetNonUnique] = movingControlPointIndices[i];
            }

            reduce( globalStaticPointsList, sumOp<labelList>() );
            reduce( globalFixedPointsList, sumOp<labelList>() );

            if ( not faceCellCenters )
            {
                reduce( globalMovingPointsList, sumOp<labelList>() );
                reduce( globalMovingPointsPatchIds, sumOp<labelList>() );
                reduce( globalMovingPointsIndices, sumOp<labelList>() );
            }

            // Construct a list of static control points which indicate whether
            // should be included or not.

            globalStaticPointsListEnabled.resize( nbStaticFaceCenters );
            globalStaticPointsListEnabled = 0;
            globalFixedPointsListEnabled.resize( nbFixedFaceCenters );
            globalFixedPointsListEnabled = 0;
            globalMovingPointsListEnabled.resize( nbMovingFaceCenters );
            globalMovingPointsListEnabled = 0;
            forAll( globalStaticPointsList, i )
            {
                // Only add the static vertex point if it's not already added to the list
                std::vector<int>::iterator it;

                it = std::find( staticControlGlobalPointLabels.begin(), staticControlGlobalPointLabels.end(), globalStaticPointsList[i] );

                if ( it == staticControlGlobalPointLabels.end() || staticControlGlobalPointLabels.size() == 0 )
                {
                    staticControlGlobalPointLabels.push_back( globalStaticPointsList[i] );
                    globalStaticPointsListEnabled[i] = 1;
                }
            }

            forAll( globalFixedPointsList, i )
            {
                // Only add the static vertex point if it's not already added to the list
                std::vector<int>::iterator it;

                it = std::find( staticControlGlobalPointLabels.begin(), staticControlGlobalPointLabels.end(), globalFixedPointsList[i] );

                if ( it == staticControlGlobalPointLabels.end() || staticControlGlobalPointLabels.size() == 0 )
                {
                    it = std::find( fixedControlGlobalPointLabels.begin(), fixedControlGlobalPointLabels.end(), globalFixedPointsList[i] );

                    if ( it == fixedControlGlobalPointLabels.end() || fixedControlGlobalPointLabels.size() == 0 )
                    {
                        fixedControlGlobalPointLabels.push_back( globalFixedPointsList[i] );
                        globalFixedPointsListEnabled[i] = 1;
                    }
                }
            }

            if ( not faceCellCenters )
            {
                forAll( movingPatchIDs, patchI )
                {
                    const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchI]].meshPoints();
                    globalMovingPointsLabelList[movingPatchIDs[patchI]] = labelList( meshPoints.size(), 0 );
                }

                forAll( globalMovingPointsList, i )
                {
                    // Only add the static vertex point if it's not already added to the list
                    std::vector<int>::iterator it;

                    it = std::find( staticControlGlobalPointLabels.begin(), staticControlGlobalPointLabels.end(), globalMovingPointsList[i] );

                    if ( it == staticControlGlobalPointLabels.end() || staticControlGlobalPointLabels.size() == 0 )
                    {
                        it = std::find( fixedControlGlobalPointLabels.begin(), fixedControlGlobalPointLabels.end(), globalMovingPointsList[i] );

                        if ( it == fixedControlGlobalPointLabels.end() || fixedControlGlobalPointLabels.size() == 0 )
                        {
                            it = std::find( movingControlGlobalPointLabels.begin(), movingControlGlobalPointLabels.end(), globalMovingPointsList[i] );

                            if ( it == movingControlGlobalPointLabels.end() || movingControlGlobalPointLabels.size() == 0 )
                            {
                                movingControlGlobalPointLabels.push_back( globalMovingPointsList[i] );
                                globalMovingPointsListEnabled[i] = 1;

                                if ( i < static_cast<int>( movingControlPointLabels.size() ) + globalMovingOffsetNonUnique
                                    && i >= globalMovingOffsetNonUnique )
                                {
                                    label patchId = globalMovingPointsPatchIds[i];
                                    label index = globalMovingPointsIndices[i];
                                    globalMovingPointsLabelList[patchId][index] = 1;
                                }
                            }
                        }
                    }
                }
            }

            // Count the number of local unique static points
            nbStaticFaceCenters = 0;

            for ( unsigned int i = 0; i < staticControlPointLabels.size(); i++ )
            {
                if ( globalStaticPointsListEnabled[i + globalStaticOffsetNonUnique] == 1 )
                    nbStaticFaceCenters++;
            }

            nbFixedFaceCenters = 0;

            for ( unsigned int i = 0; i < fixedControlPointLabels.size(); i++ )
            {
                if ( globalFixedPointsListEnabled[i + globalFixedOffsetNonUnique] == 1 )
                    nbFixedFaceCenters++;
            }

            if ( not faceCellCenters )
            {
                nbMovingFaceCenters = 0;

                for ( unsigned int i = 0; i < movingControlPointLabels.size(); i++ )
                {
                    if ( globalMovingPointsListEnabled[i + globalMovingOffsetNonUnique] == 1 )
                        nbMovingFaceCenters++;
                }
            }
        }

        nbMovingFaceCenters = 4;

        // Calculate sum of all faces on each processor
        nbGlobalStaticFaceCenters = 0;
        nbGlobalFixedFaceCenters = 0;
        nbGlobalMovingFaceCenters = 0;
        nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
        nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
        nbGlobalFixedFaceCenters[Pstream::myProcNo()] = nbFixedFaceCenters;
        nbGlobalFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters + nbStaticFaceCenters + nbFixedFaceCenters;

        reduce( nbGlobalMovingFaceCenters, sumOp<labelList>() );
        reduce( nbGlobalStaticFaceCenters, sumOp<labelList>() );
        reduce( nbGlobalFixedFaceCenters, sumOp<labelList>() );
        reduce( nbGlobalFaceCenters, sumOp<labelList>() );
    }

    nbMovingFaceCenters = sum( nbGlobalMovingFaceCenters );
    nbStaticFaceCenters = sum( nbGlobalStaticFaceCenters );
    nbFixedFaceCenters = sum( nbGlobalFixedFaceCenters );
    nbFaceCenters = sum( nbGlobalFaceCenters );

    // Determine the offset taking into account multiple processors

    int globalMovingOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalMovingOffset += nbGlobalMovingFaceCenters[i];

    int globalStaticOffset = nbMovingFaceCenters;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalStaticOffset += nbGlobalStaticFaceCenters[i];

    int globalFixedOffset = nbMovingFaceCenters + nbStaticFaceCenters;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalFixedOffset += nbGlobalFixedFaceCenters[i];

    if ( !rbf->rbf->computed )
    {
        rbf::matrix positions( nbFaceCenters, mesh().nGeometricD() );
        positions.setZero();
        int index = 0;

        const Foam::pointField & points = mesh().points();

        vectorField positionsField( positions.rows(), vector::zero );

        if ( faceCellCenters )
        {
            int offset = 0;

            forAll( movingPatchIDs, i )
            {
                const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

                // Set the positions for patch i
                forAll( faceCentres, j )
                {
                    positionsField[j + offset + globalMovingOffset] = faceCentres[j];
                }

                offset += faceCentres.size();
            }
        }

        if ( not faceCellCenters )
        {
            for ( unsigned int i = 0; i < movingControlPointLabels.size(); i++ )
            {
                if ( globalMovingPointsListEnabled[i + globalMovingOffsetNonUnique] == 1 )
                {
                    positionsField[index + globalMovingOffset] = points[movingControlPointLabels[i]];
                    index++;
                }
            }

            index = 4;
            assert( index == nbGlobalMovingFaceCenters[Pstream::myProcNo()] );
        }

        index = 0;

        for ( unsigned int i = 0; i < staticControlPointLabels.size(); i++ )
        {
            if ( globalStaticPointsListEnabled[i + globalStaticOffsetNonUnique] == 1 )
            {
                positionsField[index + globalStaticOffset] = points[staticControlPointLabels[i]];
                index++;
            }
        }

        assert( index == nbGlobalStaticFaceCenters[Pstream::myProcNo()] );

        index = 0;

        for ( unsigned int i = 0; i < fixedControlPointLabels.size(); i++ )
        {
            if ( globalFixedPointsListEnabled[i + globalFixedOffsetNonUnique] == 1 )
            {
                positionsField[index + globalFixedOffset] = points[fixedControlPointLabels[i]];
                index++;
            }
        }

        assert( index == nbGlobalFixedFaceCenters[Pstream::myProcNo()] );

        reduce( positionsField, sumOp<vectorField>() );

        // Copy the FOAM vector field to an Eigen matrix
        for ( int i = 0; i < positions.rows(); i++ )
            for ( int j = 0; j < positions.cols(); j++ )
                positions( i, j ) = positionsField[i][j];

        //std::cout << "positions = \n" << positions << std::endl;

        positions.row(0) << 0.5, 0.5;
        positions.row(1) << 0.52, 0.5;
        positions.row(2) << 0.52, 0.52;
        positions.row(3) << 0.5, 0.52;

        //std::cout << "positions = \n" << positions << std::endl;

        /*
         * Step 2: Build a matrix with the positions of every vertex in the local mesh.
         * This is only local information and does not need to be communicated to other
         * processors.
         */

        // Determine the number of points by using the 2d corrector
        nbPoints = 0;
        forAll( points, i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
                nbPoints++;
        }

        rbf::matrix positionsInterpolation( nbPoints, positions.cols() );

        index = 0;
        forAll( points, i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
            {
                for ( int j = 0; j < positionsInterpolation.cols(); j++ )
                    positionsInterpolation( index, j ) = points[i][j];

                index++;
            }
        }

        rbf->compute( positions, positionsInterpolation );

        rbf->setNbMovingAndStaticFaceCenters( nbMovingFaceCenters, nbStaticFaceCenters + nbFixedFaceCenters );
    }

    /*
     * Step 3: Build a matrix with the displacement/motion of the face center
     * positions of the static patches and the moving patches.
     * The motion needs to be communicated to every process at every mesh deformation.
     * This is considered to be the most expensive step with regards to parallel
     * scalability of the overall algorithm.
     */

    rbf::matrix values( nbFaceCenters, mesh().nGeometricD() );
    values.setZero();

    vectorField valuesField( values.rows(), vector::zero );

    if ( faceCellCenters )
    {
        int offset = 0;

        forAll( movingPatchIDs, i )
        {
            const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

            forAll( motionCenters[movingPatchIDs[i]], j )
            {
                valuesField[j + offset + globalMovingOffset] = motionCenters[movingPatchIDs[i]][j];
            }

            offset += faceCentres.size();
        }
    }

    valuesField[0] = calcVelocity();

    for ( int i = 1; i < 4; i++ )
        valuesField[i] = valuesField[0];

    if ( not faceCellCenters )
    {
        int index = 0;

        forAll( movingPatchIDs, i )
        {
            forAll( globalMovingPointsLabelList[movingPatchIDs[i]], j )
            {
                if ( globalMovingPointsLabelList[movingPatchIDs[i]][j] == 1 )
                {
                    valuesField[index + globalMovingOffset] = motionCenters[movingPatchIDs[i]][j];
                    index++;
                }
            }
        }

        index = 4;
        assert( index == nbGlobalMovingFaceCenters[Pstream::myProcNo()] );
    }

    reduce( valuesField, sumOp<vectorField>() );

    values.setZero();

    // Copy the FOAM vector field to an Eigen matrix
    for ( int i = 0; i < values.rows(); i++ )
        for ( int j = 0; j < values.cols(); j++ )
            values( i, j ) = valuesField[i][j];

    //std::cout << "values = \n" << values << std::endl;

    /*
     * Step 4: Perform the interpolation from the face centers to the complete mesh
     */

    rbf::matrix valuesInterpolation( nbPoints, values.cols() );
    valuesInterpolation.setZero();

    if ( cpu )
        rbf->rbf->computed = false;

    rbf->interpolate( values, valuesInterpolation );

    // Apply the 2d correction

    vectorField valuesInterpolationField( mesh().points().size(), Foam::vector::zero );
    int index = 0;
    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] == 0 )
        {
            for ( int j = 0; j < valuesInterpolation.cols(); j++ )
                valuesInterpolationField[i][j] = valuesInterpolation( index, j );

            index++;
        }
    }

    twoDCorrector.setShadowSide( valuesInterpolationField );

    /*
     * Step 5: Correct the mesh vertices of the fixed patches. Set these displacements to zero.
     */

    // Loop over all the patches, and set the fixed patches to zero.

    forAll( mesh().boundaryMesh(), i )
    {
        const labelList & meshPoints = mesh().boundaryMesh()[i].meshPoints();

        bool isFixedPatch = false;
        forAll( fixedPatchIDs, j )
        {
            if ( i == fixedPatchIDs[j] )
                isFixedPatch = true;
        }

        if ( isFixedPatch )
        {
            for ( int j = 0; j < meshPoints.size(); j++ )
                valuesInterpolationField[meshPoints[j]] = Foam::vector::zero;
        }
    }

    /*
     * Step 6: Set the motion of the mesh vertices
     */

    assert( newPoints.size() == valuesInterpolationField.size() );

    newPoints = valuesInterpolationField;
}
