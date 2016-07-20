
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
    Uinitial
    (
        IOobject
        (
            "U",
            runTime->timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
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
    twoDCorrector( mesh )
{
    // Find IDs of staticPatches_
    forAll( movingPatches, patchI )
    {
        label patchIndex = mesh.boundaryMesh().findPatchID( movingPatches[patchI] );

        assert( patchIndex >= 0 );

        movingPatchIDs[patchI] = patchIndex;
    }

    N = getInterfaceSize();
    dim = mesh.nGeometricD();
    data.resize( N, dim );
    data.setZero();

    readCouplingProperties();

    assert( readLabel( runTime->controlDict().lookup( "writePrecision" ) ) >= 12 );
}

foamSolidSolver::~foamSolidSolver()
{}

void foamSolidSolver::finalizeTimeStep()
{
    assert( init );

    init = false;
}

void foamSolidSolver::getDisplacementLocal( matrix & displacement )
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

    assert( movingPoints.size() == movingPointLabels.size() );

    displacement.resize( movingPoints.size(), mesh.nGeometricD() );

    for ( int i = 0; i < displacement.rows(); i++ )
        for ( int j = 0; j < displacement.cols(); j++ )
            displacement( i, j ) = pointU[movingPointLabels[i]][j];
}

int foamSolidSolver::getInterfaceSize()
{
    if ( sum( nGlobalCenters ) > 0 )
        return sum( nGlobalCenters );

    nGlobalCenters = 0;

    nGlobalCenters[Pstream::myProcNo()] = getInterfaceSizeLocal();

    reduce( nGlobalCenters, sumOp<labelList>() );

    return sum( nGlobalCenters );
}

int foamSolidSolver::getInterfaceSizeLocal()
{
    matrix pos;
    getWritePositions( pos );
    return pos.rows();
}

void foamSolidSolver::getReadPositions( matrix & readPositions )
{
    vectorField readPositionsField( getInterfaceSize(), Foam::vector::zero );

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
    getWritePositionsLocal( readPositions );
}

void foamSolidSolver::getWritePositions( matrix & writePositions )
{
    // Write positions: the face vertices of the moving patches

    movingPoints.clear();
    movingPointLabels.clear();

    int totalPoints = 0;

    forAll( movingPatchIDs, patchI )
    {
        const labelList & meshPoints = mesh.boundaryMesh()[movingPatchIDs[patchI]].meshPoints();

        totalPoints += meshPoints.size();

        forAll( meshPoints, j )
        {
            if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                continue;

            if ( movingPoints.find( meshPoints[j] ) == movingPoints.end() )
            {
                std::unordered_map<std::string, unsigned int> movingPoint;
                movingPoint["patch-id"] = movingPatchIDs[patchI];
                movingPoint["point-index"] = j;
                movingPoint["id"] = movingPoints.size() - 1;
                movingPoints[meshPoints[j]] = movingPoint;
                movingPointLabels.push_back( meshPoints[j] );
            }
        }
    }

    assert( movingPointLabels.size() == movingPoints.size() );

    writePositions.resize( movingPoints.size(), mesh.nGeometricD() );

    const Foam::pointField & points = mesh.points();

    for ( int i = 0; i < writePositions.rows(); i++ )
        for ( int j = 0; j < writePositions.cols(); j++ )
            writePositions( i, j ) = points[movingPointLabels[i]][j];
}

void foamSolidSolver::getWritePositionsLocal( matrix & writePositions )
{
    writePositions.resize( getInterfaceSizeLocal(), mesh.nGeometricD() );

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        const vectorField faceCentres( mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres() );

        for ( int i = 0; i < faceCentres.size(); i++ )
            for ( int j = 0; j < writePositions.cols(); j++ )
                writePositions( i + offset, j ) = faceCentres[i][j] + Uinitial.boundaryField()[movingPatchIDs[patchI]][i][j];

        offset += faceCentres.size();
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
    assert( traction.rows() == getInterfaceSizeLocal() );
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

    std::clock_t t = std::clock();

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

    vectorField outputField( getInterfaceSize(), Foam::vector::zero );

    matrix displacementLocal;
    getDisplacementLocal( displacementLocal );

    for ( int i = 0; i < displacementLocal.rows(); i++ )
        for ( int j = 0; j < displacementLocal.cols(); j++ )
            outputField[i + globalOffset][j] = displacementLocal( i, j );

    reduce( outputField, sumOp<vectorField>() );

    output.resize( outputField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < output.rows(); i++ )
        for ( int j = 0; j < output.cols(); j++ )
            output( i, j ) = outputField[i][j];

    assert( output.rows() == data.rows() );
    assert( output.cols() == data.cols() );

    data = output;

    if ( debug > 0 )
    {
        t = std::clock() - t;
        scalar runTime = static_cast<scalar>(t) / CLOCKS_PER_SEC;
        totalRunTime += runTime;
        totalNbIterations++;
        Info << "runtime = " << runTime << " s" << endl;
        Info << "average runtime = " << totalRunTime / totalNbIterations << " s" << endl;
        Info << "total runtime = " << totalRunTime << " s" << endl;
    }
}
