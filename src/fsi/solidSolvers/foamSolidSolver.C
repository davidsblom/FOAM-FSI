
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "foamSolidSolver.H"

Foam::debug::debugSwitch foamSolidSolver::debug( "FsiSolver", 0 );

foamSolidSolver::foamSolidSolver (
    string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    )
    :
    BaseMultiLevelSolver( 5, 3, 0 ),
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
    totalNbIterations( 0 )
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
    displacement.resize( getInterfaceSizeLocal(), mesh.nGeometricD() );

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = U.boundaryField()[movingPatchIDs[patchI]].size();

        for ( int i = 0; i < size; i++ )
            for ( int j = 0; j < displacement.cols(); j++ )
                displacement( i + offset, j ) = U.boundaryField()[movingPatchIDs[patchI]][i][j] - Uinitial.boundaryField()[movingPatchIDs[patchI]][i][j];

        offset += size;
    }
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
    int size = 0;

    forAll( movingPatchIDs, patchI )
    {
        size += mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();
    }

    return size;
}

void foamSolidSolver::getReadPositions( matrix & readPositions )
{
    getWritePositions( readPositions );
}

void foamSolidSolver::getReadPositionsLocal( matrix & readPositions )
{
    getWritePositionsLocal( readPositions );
}

void foamSolidSolver::getWritePositions( matrix & writePositions )
{
    vectorField writePositionsField( getInterfaceSize(), Foam::vector::zero );

    int globalOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalOffset += nGlobalCenters[i];

    matrix writePositionsLocal;
    getWritePositionsLocal( writePositionsLocal );

    for ( int i = 0; i < writePositionsLocal.rows(); i++ )
        for ( int j = 0; j < writePositionsLocal.cols(); j++ )
            writePositionsField[i + globalOffset][j] = writePositionsLocal( i, j );

    reduce( writePositionsField, sumOp<vectorField>() );

    writePositions.resize( writePositionsField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < writePositions.rows(); i++ )
        for ( int j = 0; j < writePositions.cols(); j++ )
            writePositions( i, j ) = writePositionsField[i][j];
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
        double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
        totalRunTime += runTime;
        totalNbIterations++;
        Info << "runtime = " << runTime << " s" << endl;
        Info << "average runtime = " << totalRunTime / totalNbIterations << " s" << endl;
        Info << "total runtime = " << totalRunTime << " s" << endl;
    }
}
