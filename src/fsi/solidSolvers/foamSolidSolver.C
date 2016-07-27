
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "foamSolidSolver.H"
#include "leastSquaresVolPointInterpolation.H"

Foam::debug::debugSwitch foamSolidSolver::debug( "FsiSolver", 0 );

foamSolidSolver::foamSolidSolver (
    string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    )
    :
    BaseMultiLevelSolver( 5, 3 ),
    name( name ),
    args( args ),
    runTime( runTime ),
    mesh
    (
        Foam::IOobject
        (
            name,
            runTime->timeName(),
            *runTime,
            Foam::IOobject::MUST_READ
        )
    ),
    U
    (
        IOobject
        (
            "U",
            runTime->timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    couplingProperties
    (
        IOobject
        (
            "couplingProperties",
            runTime->constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    movingPatches( couplingProperties.lookup( "movingSolidPatches" ) ),
    movingPatchIDs( movingPatches.size() ),
    forces(),
    nGlobalCenters( Pstream::nProcs(), 0 ),
    totalRunTime( 0 ),
    totalNbIterations( 0 ),
    twoDCorrector( mesh ),
    nbGlobalPoints( Pstream::nProcs(), 0 )
{
    // Find IDs of staticPatches_
    forAll( movingPatches, patchI )
    {
        label patchIndex = mesh.boundaryMesh().findPatchID( movingPatches[patchI] );

        assert( patchIndex >= 0 );

        movingPatchIDs[patchI] = patchIndex;
    }

    matrix writePositions;
    getWritePositions( writePositions );

    N = writePositions.rows();

    dim = mesh.nGeometricD();
    data.resize( N, dim );
    data.setZero();

    readCouplingProperties();

    assert( readLabel( runTime->controlDict().lookup( "writePrecision" ) ) >= 12 );
    timeIndex = runTime->timeIndex();
}

foamSolidSolver::~foamSolidSolver()
{}

void foamSolidSolver::finalizeTimeStep()
{
    assert( init );

    init = false;
}

void foamSolidSolver::getDisplacementLocal( matrix & displacementLocal )
{
    leastSquaresVolPointInterpolation pointInterpolation( mesh );

    // Create point mesh
    pointMesh pMesh( mesh );

    wordList types
    (
        pMesh.boundary().size(),
        calculatedFvPatchVectorField::typeName
    );

    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime->timeName(),
            mesh
        ),
        pMesh,
        dimensionedVector( "zero", dimLength, Foam::vector::zero ),
        types
    );

    pointInterpolation.interpolate( U, pointU );

    assert( globalMovingPoints.size() == globalMovingPointLabels.size() );

    unsigned int N = 0;

    for ( unsigned i = 0; i < globalMovingPoints.size(); i++ )
        if ( globalMovingPoints[globalMovingPointLabels[i]]["local-point"] == true )
            N++;

    displacementLocal.resize( N, mesh.nGeometricD() );

    unsigned int i = 0;

    for ( auto point : globalMovingPoints )
    {
        bool localPoint = point.second["local-point"];

        if ( localPoint )
        {
            int localId = point.second["local-id"];

            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                displacementLocal( i, j ) = pointU[movingPointLabels[localId]][j];
            }

            i++;
        }
    }
}

void foamSolidSolver::getReadPositions( matrix & readPositions )
{
    int size = 0;

    forAll( movingPatchIDs, patchI )
    {
        size += mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();
    }

    nGlobalCenters = 0;
    nGlobalCenters[Pstream::myProcNo()] = size;
    reduce( nGlobalCenters, sumOp<labelList>() );

    vectorField readPositionsField( sum( nGlobalCenters ), Foam::vector::zero );

    int globalOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalOffset += nGlobalCenters[i];

    matrix readPositionsLocal;
    getReadPositionsLocal( readPositionsLocal );

    for ( int i = 0; i < readPositionsLocal.rows(); i++ )
        for ( int j = 0; j < readPositionsLocal.cols(); j++ )
            readPositionsField[i + globalOffset][j] = readPositionsLocal( i, j );

    reduce( readPositionsField, sumOp<vectorField>() );

    readPositions.resize( readPositionsField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < readPositions.rows(); i++ )
        for ( int j = 0; j < readPositions.cols(); j++ )
            readPositions( i, j ) = readPositionsField[i][j];
}

void foamSolidSolver::getReadPositionsLocal( matrix & readPositions )
{
    int size = 0;

    forAll( movingPatchIDs, patchId )
    {
        size += mesh.boundaryMesh()[movingPatchIDs[patchId]].faceCentres().size();
    }

    readPositions.resize( size, mesh.nGeometricD() );

    int offset = 0;

    forAll( movingPatchIDs, patchId )
    {
        const vectorField faceCentres( mesh.boundaryMesh()[movingPatchIDs[patchId]].faceCentres() );

        for ( int i = 0; i < faceCentres.size(); i++ )
            for ( int j = 0; j < readPositions.cols(); j++ )
                readPositions( i + offset, j ) = faceCentres[i][j];

        offset += faceCentres.size();
    }
}

void foamSolidSolver::getWritePositions( matrix & writePositions )
{
    // Write positions: the face vertices of the moving patches

    movingPoints.clear();
    movingPointLabels.clear();

    globalMovingPoints.clear();
    globalMovingPointLabels.clear();

    forAll( movingPatchIDs, patchI )
    {
        const labelList & meshPoints = mesh.boundaryMesh()[movingPatchIDs[patchI]].meshPoints();

        forAll( meshPoints, j )
        {
            if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                continue;

            if ( movingPoints.find( meshPoints[j] ) == movingPoints.end() )
            {
                std::unordered_map<std::string, unsigned int> movingPoint;
                movingPoint["patch-id"] = movingPatchIDs[patchI];
                movingPoint["point-index"] = j;
                movingPoint["local-id"] = movingPoints.size();
                movingPoint["mesh-index"] = meshPoints[j];
                movingPoints[meshPoints[j]] = movingPoint;
                movingPointLabels.push_back( meshPoints[j] );
            }
        }
    }

    assert( movingPointLabels.size() == movingPoints.size() );

    if ( Pstream::nProcs() > 1 )
    {
        IOobject addrHeader
        (
            "pointProcAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        );

        assert( addrHeader.headerOk() );
        labelIOList pointProcAddressing( addrHeader );

        assert( pointProcAddressing.size() == mesh.points().size() );

        nbGlobalPoints = 0;
        nbGlobalPoints[Pstream::myProcNo()] = movingPointLabels.size();
        reduce( nbGlobalPoints, sumOp<labelList>() );

        // Construct a list with all the global point labels
        // Thereafter, construct a list which indicates whether the point is
        // already included or not. Use this later to build a list of the
        // unique points.
        labelList globalPointsIndices( sum( nbGlobalPoints ), 0 );

        int globalOffsetNonUnique = 0;

        for ( int i = 0; i < Pstream::myProcNo(); i++ )
            globalOffsetNonUnique += nbGlobalPoints[i];

        for ( auto label : movingPoints )
        {
            globalPointsIndices[label.second["local-id"] + globalOffsetNonUnique] = pointProcAddressing[label.first];
        }

        reduce( globalPointsIndices, sumOp<labelList>() );

        // Construct unique list with global indices

        forAll( globalPointsIndices, i )
        {
            if ( globalMovingPoints.find( globalPointsIndices[i] ) == globalMovingPoints.end() )
            {
                std::unordered_map<std::string, unsigned int> movingPoint;
                movingPoint["global-id"] = globalMovingPoints.size();
                int localId = i - globalOffsetNonUnique;

                if ( localId >= 0
                    && !movingPoints.empty()
                    && localId < static_cast<int>( movingPoints.size() ) )
                {
                    movingPoint["local-point"] = true;
                    movingPoint["mesh-index"] = movingPointLabels[localId];
                    movingPoint["local-id"] = localId;
                    movingPoint["patch-id"] = movingPoints[movingPoint["mesh-index"]]["patch-id"];
                    movingPoint["point-index"] = movingPoints[movingPoint["mesh-index"]]["point-index"];
                }
                else
                {
                    movingPoint["local-point"] = false;
                }

                globalMovingPoints[globalPointsIndices[i]] = movingPoint;
                globalMovingPointLabels.push_back( globalPointsIndices[i] );
            }
        }
    }
    else
    {
        for ( auto point : movingPoints )
        {
            point.second["local-point"] = true;
            point.second["global-id"] = point.second["local-id"];
            globalMovingPoints[point.first] = point.second;
        }

        for ( auto label : movingPointLabels )
        {
            globalMovingPointLabels.push_back( label );
        }
    }

    vectorField positionsField( globalMovingPoints.size(), Foam::vector::zero );

    for ( auto point : globalMovingPoints )
    {
        if ( point.second["local-point"] )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                positionsField[point.second["global-id"]][j] = mesh.points()[point.second["mesh-index"]][j];
            }
        }
    }

    reduce( positionsField, sumOp<vectorField>() );

    writePositions.resize( positionsField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < writePositions.rows(); i++ )
        for ( int j = 0; j < writePositions.cols(); j++ )
            writePositions( i, j ) = positionsField[i][j];
}

void foamSolidSolver::getWritePositionsLocal( matrix & writePositionsLocal )
{
    matrix writePositions;
    getWritePositions( writePositions );

    unsigned int N = 0;

    for ( auto point : globalMovingPoints )
        if ( point.second["local-point"] == true )
            N++;

    writePositionsLocal.resize( N, mesh.nGeometricD() );

    int i = 0;

    for ( auto point : globalMovingPoints )
    {
        if ( point.second["local-point"] == true )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                writePositionsLocal( i, j ) = mesh.points()[point.second["mesh-index"]][j];
            }

            i++;
        }
    }
}

void foamSolidSolver::readCouplingProperties()
{
    // Grab solid patch field

    forAll( movingPatchIDs, patchI )
    {
        assert( U.boundaryField()[movingPatchIDs[patchI]].type() == solidTractionFvPatchVectorField::typeName );

        solidTractionFvPatchVectorField & tForce =
            refCast<solidTractionFvPatchVectorField>
            (
            U.boundaryField()[movingPatchIDs[patchI]]
            );

        forces.push_back( &tForce );
    }
}

void foamSolidSolver::setTractionLocal( const matrix & traction )
{
    assert( traction.cols() == mesh.nGeometricD() );

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        vectorField tractionSolid( size, Foam::vector::zero );

        for ( int i = 0; i < size; i++ )
            for ( int j = 0; j < traction.cols(); j++ )
                tractionSolid[i][j] = traction( i + offset, j );

        forces.at( patchI )->traction() = tractionSolid;

        offset += size;
    }
}

void foamSolidSolver::run()
{
    while ( isRunning() )
    {
        Info << "Time = " << runTime->timeName() << endl;

        initTimeStep();
        solve();
        finalizeTimeStep();
    }
}

void foamSolidSolver::solve(
    const matrix & input,
    matrix & output
    )
{
    assert( input.rows() == sum( nGlobalCenters ) );
    assert( input.cols() == mesh.nGeometricD() );

    int globalOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalOffset += nGlobalCenters[i];

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        vectorField tractionSolid( size, Foam::vector::zero );

        for ( int i = 0; i < size; i++ )
            for ( int j = 0; j < input.cols(); j++ )
                tractionSolid[i][j] = input( i + offset + globalOffset, j );

        forces.at( patchI )->traction() = tractionSolid;

        offset += size;
    }

    solve();

    leastSquaresVolPointInterpolation pointInterpolation( mesh );

    // Create point mesh
    pointMesh pMesh( mesh );

    wordList types
    (
        pMesh.boundary().size(),
        calculatedFvPatchVectorField::typeName
    );

    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime->timeName(),
            mesh
        ),
        pMesh,
        dimensionedVector( "zero", dimLength, Foam::vector::zero ),
        types
    );

    pointInterpolation.interpolate( U, pointU );

    assert( globalMovingPoints.size() == globalMovingPointLabels.size() );

    vectorField outputField( globalMovingPoints.size(), Foam::vector::zero );

    for ( unsigned i = 0; i < globalMovingPoints.size(); i++ )
    {
        bool localPoint = globalMovingPoints[globalMovingPointLabels[i]]["local-point"];

        if ( localPoint )
        {
            int localId = globalMovingPoints[globalMovingPointLabels[i]]["local-id"];

            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                outputField[i][j] = pointU[movingPointLabels[localId]][j];
            }
        }
    }

    reduce( outputField, sumOp<vectorField>() );

    output.resize( outputField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < output.rows(); i++ )
        for ( int j = 0; j < output.cols(); j++ )
            output( i, j ) = outputField[i][j];

    data = output;
}
