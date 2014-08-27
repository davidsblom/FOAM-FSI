
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "foamSolidSolver.H"

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
  solidPatchName( couplingProperties.lookup( "solidPatch" ) ),
  solidPatchID( mesh.boundaryMesh().findPatchID( solidPatchName ) ),
  nPoints( Pstream::nProcs(), 0 ),
  nGlobalPoints( Pstream::nProcs(), 0 ),
  nGlobalCenters( Pstream::nProcs(), 0 ),
  globalPointsUnique( 1, -1 ),
  globalPointsNonUnique( 1, -1 )
{
  assert( solidPatchID >= 0 );

  N = getInterfacePointsSize();
  dim = 3;
  data.resize( N, dim );
  data.setZero();

  readCouplingProperties();
}

foamSolidSolver::~foamSolidSolver()
{}

void foamSolidSolver::finalizeTimeStep()
{
  assert( init );

  if ( name == fvMesh::defaultRegion )
  {
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;
  }

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

  pointVectorField pointDU
  (
    IOobject
    (
      "pointDU",
      runTime->timeName(),
      mesh
    ),
    pMesh,
    dimensionedVector( "zero", dimLength, Foam::vector::zero ),
    types
  );

  pointInterpolation.interpolate( U, pointDU );

  const labelList & meshPoints = mesh.boundaryMesh()[solidPatchID].meshPoints();

  vectorField outputFieldNonUnique( globalPointsNonUnique.size(), Foam::vector::zero );
  vectorField outputFieldUnique( globalPointsUnique.size(), Foam::vector::zero );

  assert( sum( nPoints ) == globalPointsUnique.size() );
  assert( sum( nGlobalPoints ) == globalPointsNonUnique.size() );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalPoints[i];

  forAll( meshPoints, i )
  {
    outputFieldNonUnique[i + offset] = pointDU[meshPoints[i]];
  }

  reduce( outputFieldNonUnique, sumOp<vectorField>() );

  forAll( outputFieldNonUnique, i )
  {
    outputFieldUnique[globalPointsNonUnique[i]] = outputFieldNonUnique[i];
  }

  displacement.resize( nPoints[Pstream::myProcNo()], 3 );

  offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nPoints[i];

  for ( int i = 0; i < displacement.rows(); i++ )
    for ( int j = 0; j < displacement.cols(); j++ )
      displacement( i, j ) = outputFieldUnique[i + offset][j];
}

int foamSolidSolver::getInterfaceCentersSize()
{
  if ( sum( nGlobalCenters ) > 0 )
    return sum( nGlobalCenters );

  nGlobalCenters = 0;
  nGlobalCenters[Pstream::myProcNo()] = mesh.boundaryMesh()[solidPatchID].size();

  reduce( nGlobalCenters, sumOp<labelList>() );

  return sum( nGlobalCenters );
}

int foamSolidSolver::getInterfaceCentersSizeLocal()
{
  return mesh.boundaryMesh()[solidPatchID].size();
}

int foamSolidSolver::getInterfacePointsSize()
{
  if ( sum( nPoints ) > 0 )
    return sum( nPoints );

  const pointField & points = mesh.points();
  labelList markedPoints( points.size(), 0 );

  const polyBoundaryMesh & bMesh( mesh.boundaryMesh() );
  const labelList & procPatches( mesh.globalData().processorPatches() );

  forAll( procPatches, i )
  {
    const processorPolyPatch & procPatch = refCast<const processorPolyPatch>( bMesh[procPatches[i]] );
    const labelList & procPatchPoints = procPatch.meshPoints();

    if ( procPatch.neighbour() )    // If not the owner of the patch
    {
      forAll( procPatchPoints, k )    // Go over all the points and mark them with -5
      {
        markedPoints[procPatchPoints[k]] = -5;
      }
    }
  }

  // Count the owned points
  label nbPoints = 0;

  const labelList & meshPoints = mesh.boundaryMesh()[solidPatchID].meshPoints();

  forAll( meshPoints, i )
  {
    if ( markedPoints[meshPoints[i]] != -5 )
      nbPoints++;
  }

  nPoints = 0;
  nPoints[Pstream::myProcNo()] = nbPoints;

  reduce( nPoints, sumOp<labelList>() );

  return sum( nPoints );
}

void foamSolidSolver::getReadPositions( matrix & readPositions )
{
  vectorField readPositionsField( getInterfaceCentersSize(), Foam::vector::zero );

  const vectorField faceCentres( mesh.boundaryMesh()[solidPatchID].faceCentres() );

  assert( faceCentres.size() == nGlobalCenters[Pstream::myProcNo()] );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalCenters[i];

  forAll( faceCentres, i )
  {
    readPositionsField[i + offset] = faceCentres[i];
  }

  reduce( readPositionsField, sumOp<vectorField>() );

  readPositions.resize( readPositionsField.size(), 3 );

  for ( int i = 0; i < readPositions.rows(); i++ )
    for ( int j = 0; j < readPositions.cols(); j++ )
      readPositions( i, j ) = readPositionsField[i][j];
}

void foamSolidSolver::getReadPositionsLocal( matrix & readPositions )
{
  const vectorField faceCentres( mesh.boundaryMesh()[solidPatchID].faceCentres() );

  readPositions.resize( faceCentres.size(), 3 );

  for ( int i = 0; i < readPositions.rows(); i++ )
    for ( int j = 0; j < readPositions.cols(); j++ )
      readPositions( i, j ) = faceCentres[i][j];
}

void foamSolidSolver::getWritePositions( matrix & writePositions )
{
  getInterfacePointsSize();

  const labelList & meshPoints = mesh.boundaryMesh()[solidPatchID].meshPoints();
  const pointField & points = mesh.points();

  // nGlobalPoints is a list containing the number of boundary points
  // on each processor core. The cores do contain double points.
  nGlobalPoints = 0;
  nGlobalPoints[Pstream::myProcNo()] = meshPoints.size();
  reduce( nGlobalPoints, sumOp<labelList>() );

  labelList globalBoundaryLabels( sum( nGlobalPoints ), 0 );

  if ( Pstream::parRun() )
  {
    // Read the processor point addressing from disk to find the
    // points which are shared by 2 or more processors.
    // This file only exists when the computation is run in
    // parallel.
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
    assert( pointProcAddressing.size() == points.size() );

    // globalBoundaryLabels contains the global mesh point labels for the
    // boundary patch.

    int offset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
      offset += nGlobalPoints[i];

    forAll( meshPoints, i )
    {
      globalBoundaryLabels[i + offset] = pointProcAddressing[meshPoints[i]];
    }

    reduce( globalBoundaryLabels, sumOp<labelList>() );
  }

  if ( !Pstream::parRun() )
  {
    forAll( meshPoints, i )
    {
      globalBoundaryLabels[i] = meshPoints[i];
    }
  }

  // Generate a list of labels to convert the unique list to the
  // double list and vice-versa.
  // The first label in globalBoundaryLabels is used for the post processing.
  // The others are copied from the first one.
  // Thus globalPointsNonUnique references the globalBoundaryLabels list.

  globalPointsNonUnique.resize( globalBoundaryLabels.size(), -1 );
  globalPointsNonUnique = -1;
  label nbUniquePoints = 0;

  forAll( globalBoundaryLabels, i )
  {
    if ( globalPointsNonUnique[i] == -1 )
    {
      // Save the index and search for doubles
      globalPointsNonUnique[i] = nbUniquePoints;

      // Search for the index in globalBoundaryLabels
      forAll( globalBoundaryLabels, j )
      {
        if ( i != j && globalBoundaryLabels[i] == globalBoundaryLabels[j] )
          globalPointsNonUnique[j] = nbUniquePoints;
      }

      nbUniquePoints++;
    }
  }

  assert( nbUniquePoints > 0 );
  assert( nbUniquePoints == sum( nPoints ) );
  forAll( globalPointsNonUnique, i )
  {
    assert( globalPointsNonUnique[i] <= nbUniquePoints );
    assert( globalPointsNonUnique[i] != -1 );
  }

  // Generate a list of labels to conver the input matrix
  // to the displacement patch with double points.

  globalPointsUnique.resize( nbUniquePoints, -1 );
  globalPointsUnique = -1;

  forAll( globalPointsNonUnique, i )
  {
    if ( globalPointsUnique[globalPointsNonUnique[i]] == -1 )
      globalPointsUnique[globalPointsNonUnique[i]] = i;
  }

  forAll( globalPointsUnique, i )
  {
    assert( globalPointsUnique[i] != -1 );
  }

  // Get the read positions from the points field

  vectorField writePositionsFieldNonUnique( globalBoundaryLabels.size(), Foam::vector::zero );
  vectorField writePositionsFieldUnique( nbUniquePoints, Foam::vector::zero );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalPoints[i];

  forAll( meshPoints, i )
  {
    writePositionsFieldNonUnique[i + offset] = points[meshPoints[i]];
  }

  reduce( writePositionsFieldNonUnique, sumOp<vectorField>() );

  forAll( writePositionsFieldNonUnique, i )
  {
    writePositionsFieldUnique[globalPointsNonUnique[i]] = writePositionsFieldNonUnique[i];
  }

  // Copy the read positions from the vector field to the Eigen matrix

  writePositions.resize( writePositionsFieldUnique.size(), 3 );

  for ( int i = 0; i < writePositions.rows(); i++ )
    for ( int j = 0; j < writePositions.cols(); j++ )
      writePositions( i, j ) = writePositionsFieldUnique[i][j];
}

void foamSolidSolver::getWritePositionsLocal( matrix & writePositions )
{
  matrix globalWritePositions;
  getWritePositions( globalWritePositions );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nPoints[i];

  writePositions.resize( nPoints[Pstream::myProcNo()], globalWritePositions.cols() );

  for ( int i = 0; i < writePositions.rows(); i++ )
    for ( int j = 0; j < writePositions.cols(); j++ )
      writePositions( i, j ) = globalWritePositions( i + offset, j );
}

void foamSolidSolver::readCouplingProperties()
{
  // Grab solid patch field
  if
  (
    U.boundaryField()[solidPatchID].type()
    != solidTractionFvPatchVectorField::typeName
  )
  {
    FatalErrorIn( args->executable() )
    << "Boundary condition on " << U.name()
    << " is "
    << U.boundaryField()[solidPatchID].type()
    << "for fluid -solid interface patch, instead "
    << solidTractionFvPatchVectorField::typeName
    << abort( FatalError );
  }

  solidTractionFvPatchVectorField & tForce =
    refCast<solidTractionFvPatchVectorField>
    (
    U.boundaryField()[solidPatchID]
    );

  this->tForce = &tForce;
}

void foamSolidSolver::setTractionLocal( const matrix & traction )
{
  assert( traction.rows() == mesh.boundaryMesh()[solidPatchID].size() );
  assert( traction.cols() == 3 );

  vectorField tractionSolid( traction.rows(), Foam::vector::zero );

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      tractionSolid[i][j] = traction( i, j );

  tForce->traction() = tractionSolid;
}

void foamSolidSolver::solve(
  const matrix & input,
  matrix & output
  )
{
  assert( input.rows() == sum( nGlobalCenters ) );
  assert( input.cols() == 3 );

  vectorField readPositionsField( getInterfaceCentersSize(), Foam::vector::zero );

  const vectorField faceCentres( mesh.boundaryMesh()[solidPatchID].faceCentres() );

  assert( faceCentres.size() == nGlobalCenters[Pstream::myProcNo()] );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalCenters[i];

  vectorField tractionSolid( mesh.boundaryMesh()[solidPatchID].size(), Foam::vector::zero );

  for ( int i = 0; i < tractionSolid.size(); i++ )
    for ( int j = 0; j < input.cols(); j++ )
      tractionSolid[i][j] = input( i + offset, j );

  tForce->traction() = tractionSolid;

  solve();

  leastSquaresVolPointInterpolation pointInterpolation( mesh );

  // Create point mesh
  pointMesh pMesh( mesh );

  wordList types
  (
    pMesh.boundary().size(),
    calculatedFvPatchVectorField::typeName
  );

  pointVectorField pointDU
  (
    IOobject
    (
      "pointDU",
      runTime->timeName(),
      mesh
    ),
    pMesh,
    dimensionedVector( "zero", dimLength, Foam::vector::zero ),
    types
  );

  pointInterpolation.interpolate( U, pointDU );

  const labelList & meshPoints = mesh.boundaryMesh()[solidPatchID].meshPoints();

  vectorField outputFieldNonUnique( globalPointsNonUnique.size(), Foam::vector::zero );
  vectorField outputFieldUnique( globalPointsUnique.size(), Foam::vector::zero );

  assert( sum( nPoints ) == globalPointsUnique.size() );
  assert( sum( nGlobalPoints ) == globalPointsNonUnique.size() );

  offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalPoints[i];

  forAll( meshPoints, i )
  {
    outputFieldNonUnique[i + offset] = pointDU[meshPoints[i]];
  }

  reduce( outputFieldNonUnique, sumOp<vectorField>() );

  forAll( outputFieldNonUnique, i )
  {
    outputFieldUnique[globalPointsNonUnique[i]] = outputFieldNonUnique[i];
  }

  output.resize( outputFieldUnique.size(), 3 );

  for ( int i = 0; i < output.rows(); i++ )
    for ( int j = 0; j < output.cols(); j++ )
      output( i, j ) = outputFieldUnique[i][j];

  assert( output.rows() == data.rows() );
  assert( output.cols() == data.cols() );

  data = output;
}
