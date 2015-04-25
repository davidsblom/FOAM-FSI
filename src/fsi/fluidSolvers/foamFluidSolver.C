
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "foamFluidSolver.H"

foamFluidSolver::foamFluidSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    )
    :
    BaseMultiLevelSolver( 5, 3, 0 ),
    name( name ),
    args( args ),
    runTime( runTime ),
    meshPtr
    (
    dynamicFvMesh::New
    (
        IOobject
        (
            name,
            runTime->timeName(),
            *runTime,
            IOobject::MUST_READ
        )
    )
    ),
    mesh( meshPtr() ),
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
    movingPatches( couplingProperties.lookup( "movingFluidPatches" ) ),
    movingPatchIDs( movingPatches.size() ),
    movingPatchesDispl( mesh.boundaryMesh().size(), vectorField( 0 ) ),
    movingPatchesDisplOld( mesh.boundaryMesh().size(), vectorField( 0 ) ),
    nPoints( Pstream::nProcs(), 0 ),
    nGlobalPoints( Pstream::nProcs(), 0 ),
    nGlobalCenters( Pstream::nProcs(), 0 ),
    globalPointsUnique( 1, -1 ),
    globalPointsNonUnique( 1, -1 ),
    totalRunTime( 0 ),
    totalNbIterations( 0 )
{
    // Find IDs of staticPatches_
    forAll( movingPatches, patchI )
    {
        label patchIndex = mesh.boundaryMesh().findPatchID( movingPatches[patchI] );

        assert( patchIndex >= 0 );

        movingPatchIDs[patchI] = patchIndex;
    }

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        movingPatchesDispl[movingPatchIDs[patchI]] = vectorField( size, Foam::vector::zero );
        movingPatchesDisplOld[movingPatchIDs[patchI]] = vectorField( size, Foam::vector::zero );
    }

    N = getInterfaceSize();
    dim = mesh.nGeometricD();
    data.resize( N, dim );
    data.setZero();
}

foamFluidSolver::~foamFluidSolver()
{}

void foamFluidSolver::finalizeTimeStep()
{
    assert( init );

    init = false;
}

int foamFluidSolver::getInterfaceSize()
{
    if ( sum( nGlobalCenters ) > 0 )
        return sum( nGlobalCenters );

    nGlobalCenters = 0;

    nGlobalCenters[Pstream::myProcNo()] = getInterfaceSizeLocal();

    reduce( nGlobalCenters, sumOp<labelList>() );

    return sum( nGlobalCenters );
}

int foamFluidSolver::getInterfaceSizeLocal()
{
    int size = 0;

    forAll( movingPatchIDs, patchI )
    {
        size += mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();
    }

    return size;
}

void foamFluidSolver::getReadPositions( matrix & readPositions )
{
    getWritePositions( readPositions );
}

void foamFluidSolver::getReadPositionsLocal( matrix & readPositions )
{
    getWritePositionsLocal( readPositions );
}

void foamFluidSolver::getWritePositions( matrix & writePositions )
{
    vectorField writePositionsField( getInterfaceSize(), Foam::vector::zero );

    int globalOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalOffset += nGlobalCenters[i];

    matrix writePositionsLocal;
    getWritePositionsLocal ( writePositionsLocal );

    for ( int i = 0; i < writePositionsLocal.rows(); i++ )
        for ( int j = 0; j < writePositionsLocal.cols(); j++ )
            writePositionsField[i + globalOffset][j] = writePositionsLocal( i, j );

    reduce( writePositionsField, sumOp<vectorField>() );

    writePositions.resize( writePositionsField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < writePositions.rows(); i++ )
        for ( int j = 0; j < writePositions.cols(); j++ )
            writePositions( i, j ) = writePositionsField[i][j];
}

void foamFluidSolver::getWritePositionsLocal( matrix & writePositions )
{
    writePositions.resize( getInterfaceSizeLocal(), mesh.nGeometricD() );

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        const vectorField faceCentres( mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres() );

        for ( int i = 0; i < faceCentres.size(); i++ )
            for ( int j = 0; j < writePositions.cols(); j++ )
                writePositions( i + offset, j ) = faceCentres[i][j];

        offset += faceCentres.size();
    }
}

void foamFluidSolver::moveMesh( Field<vectorField> motion )
{
    std::clock_t t = std::clock();

    RBFMeshMotionSolver & motionSolver =
        const_cast<RBFMeshMotionSolver &>
        (
        mesh.lookupObject<RBFMeshMotionSolver>( "dynamicMeshDict" )
        );

    motionSolver.setMotion( motion );

    mesh.update();

    t = std::clock() - t;
    double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
    Info << "timing mesh deformation = " << runTime << " s" << endl;
}

void foamFluidSolver::setDisplacementLocal( const matrix & displacement )
{
    movingPatchesDisplOld = movingPatchesDispl;

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        for ( int i = 0; i < size; i++ )
            for ( int j = 0; j < displacement.cols(); j++ )
                movingPatchesDispl[movingPatchIDs[patchI]][i][j] = displacement( i + offset, j );

        offset += size;
    }

    moveMesh( movingPatchesDisplOld - movingPatchesDispl );
}

void foamFluidSolver::run()
{
    while ( isRunning() )
    {
        initTimeStep();
        solve();
        finalizeTimeStep();
    }
}

void foamFluidSolver::solve(
    const matrix & input,
    matrix & output
    )
{
    std::clock_t t = std::clock();

    movingPatchesDisplOld = movingPatchesDispl;

    int globalOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalOffset += nGlobalCenters[i];

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        for ( int i = 0; i < size; i++ )
            for ( int j = 0; j < input.cols(); j++ )
                movingPatchesDispl[movingPatchIDs[patchI]][i][j] = input( i + offset + globalOffset, j );

        offset += size;
    }

    assert( offset == getInterfaceSizeLocal() );

    moveMesh( movingPatchesDispl - movingPatchesDisplOld );

    solve();

    matrix traction;
    getTractionLocal( traction );

    assert( traction.rows() == nGlobalCenters[Pstream::myProcNo()] );
    assert( traction.cols() == mesh.nGeometricD() );

    vectorField writePositionsField( sum( nGlobalCenters ), Foam::vector::zero );

    assert( writePositionsField.size() == data.rows() );

    for ( int i = 0; i < traction.rows(); i++ )
        for ( int j = 0; j < traction.cols(); j++ )
            writePositionsField[i + globalOffset][j] = traction( i, j );

    reduce( writePositionsField, sumOp<vectorField>() );

    output.resize( writePositionsField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < output.rows(); i++ )
        for ( int j = 0; j < output.cols(); j++ )
            output( i, j ) = writePositionsField[i][j];

    assert( output.rows() == data.rows() );
    assert( output.cols() == data.cols() );

    data = output;

    t = std::clock() - t;
    double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
    totalRunTime += runTime;
    totalNbIterations++;
    Info << "runtime = " << runTime << " s" << endl;
    Info << "average runtime = " << totalRunTime / totalNbIterations << " s" << endl;
    Info << "total runtime = " << totalRunTime << " s" << endl;
}
