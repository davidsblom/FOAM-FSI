
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
  fluidPatchName( couplingProperties.lookup( "fluidPatch" ) ),
  fluidPatchID( mesh.boundaryMesh().findPatchID( fluidPatchName ) ),
  fluidPatchPointsDispl( mesh.boundaryMesh()[fluidPatchID].meshPoints().size(), Foam::vector::zero ),
  fluidPatchPointsDisplOld( mesh.boundaryMesh()[fluidPatchID].meshPoints().size(), Foam::vector::zero ),
  nPoints( Pstream::nProcs(), 0 ),
  nGlobalPoints( Pstream::nProcs(), 0 ),
  nGlobalCenters( Pstream::nProcs(), 0 ),
  globalPointsUnique( 1, -1 ),
  globalPointsNonUnique( 1, -1 )
{
  assert( fluidPatchID >= 0 );

  N = getInterfaceCentersSize();
  dim = 3;
  data.resize( N, dim );
  data.setZero();
}

foamFluidSolver::~foamFluidSolver()
{}

void foamFluidSolver::finalizeTimeStep()
{
  assert( init );

  if ( name == fvMesh::defaultRegion || name == "fluid-level-0" )
  {
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;
  }

  init = false;
}

int foamFluidSolver::getInterfaceCentersSize()
{
  if ( sum( nGlobalCenters ) > 0 )
    return sum( nGlobalCenters );

  nGlobalCenters = 0;
  nGlobalCenters[Pstream::myProcNo()] = mesh.boundaryMesh()[fluidPatchID].size();

  reduce( nGlobalCenters, sumOp<labelList>() );

  return sum( nGlobalCenters );
}

int foamFluidSolver::getInterfacePointsSize()
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

    if ( procPatch.neighbour() )  // If not the owner of the patch
    {
      forAll( procPatchPoints, k )  // Go over all the points and mark them with -5
      {
        markedPoints[procPatchPoints[k]] = -5;
      }
    }
  }

  // Count the owned points
  label nbPoints = 0;

  const labelList & meshPoints = mesh.boundaryMesh()[fluidPatchID].meshPoints();

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

int foamFluidSolver::getInterfacePointsSizeLocal()
{
  getInterfacePointsSize();

  return nPoints[Pstream::myProcNo()];
}

void foamFluidSolver::getReadPositions( matrix & readPositions )
{
  getInterfacePointsSize();

  const labelList & meshPoints = mesh.boundaryMesh()[fluidPatchID].meshPoints();
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

  vectorField readPositionsFieldNonUnique( globalBoundaryLabels.size(), Foam::vector::zero );
  vectorField readPositionsFieldUnique( nbUniquePoints, Foam::vector::zero );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalPoints[i];

  forAll( meshPoints, i )
  {
    readPositionsFieldNonUnique[i + offset] = points[meshPoints[i]];
  }

  reduce( readPositionsFieldNonUnique, sumOp<vectorField>() );

  forAll( readPositionsFieldNonUnique, i )
  {
    readPositionsFieldUnique[globalPointsNonUnique[i]] = readPositionsFieldNonUnique[i];
  }

  // Copy the read positions from the vector field to the Eigen matrix

  readPositions.resize( readPositionsFieldUnique.size(), 3 );

  for ( int i = 0; i < readPositions.rows(); i++ )
    for ( int j = 0; j < readPositions.cols(); j++ )
      readPositions( i, j ) = readPositionsFieldUnique[i][j];
}

void foamFluidSolver::getReadPositionsLocal( matrix & readPositions )
{
  matrix globalReadPositions;
  getReadPositions( globalReadPositions );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nPoints[i];

  readPositions.resize( nPoints[Pstream::myProcNo()], globalReadPositions.cols() );

  for ( int i = 0; i < readPositions.rows(); i++ )
    for ( int j = 0; j < readPositions.cols(); j++ )
      readPositions( i, j ) = globalReadPositions( i + offset, j );
}

void foamFluidSolver::getWritePositions( matrix & writePositions )
{
  vectorField writePositionsField( getInterfaceCentersSize(), Foam::vector::zero );

  const vectorField faceCentres( mesh.boundaryMesh()[fluidPatchID].faceCentres() );

  assert( faceCentres.size() == nGlobalCenters[Pstream::myProcNo()] );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalCenters[i];

  forAll( faceCentres, i )
  {
    writePositionsField[i + offset] = faceCentres[i];
  }

  reduce( writePositionsField, sumOp<vectorField>() );

  writePositions.resize( writePositionsField.size(), 3 );

  for ( int i = 0; i < writePositions.rows(); i++ )
    for ( int j = 0; j < writePositions.cols(); j++ )
      writePositions( i, j ) = writePositionsField[i][j];
}

void foamFluidSolver::getWritePositionsLocal( matrix & writePositions )
{
  const vectorField faceCentres( mesh.boundaryMesh()[fluidPatchID].faceCentres() );

  writePositions.resize( faceCentres.size(), 3 );

  for ( int i = 0; i < writePositions.rows(); i++ )
    for ( int j = 0; j < writePositions.cols(); j++ )
      writePositions( i, j ) = faceCentres[i][j];
}

void foamFluidSolver::setDisplacementLocal( const matrix & displacement )
{
  fluidPatchPointsDisplOld = fluidPatchPointsDispl;

  vectorField inputNonUnique( globalPointsNonUnique.size(), Foam::vector::zero );
  vectorField inputUnique( globalPointsUnique.size(), Foam::vector::zero );

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nPoints[i];

  for ( int i = 0; i < displacement.rows(); i++ )
    for ( int j = 0; j < displacement.cols(); j++ )
      inputUnique[i + offset][j] = displacement( i, j );

  reduce( inputUnique, sumOp<vectorField>() );

  forAll( inputNonUnique, i )
  {
    inputNonUnique[i] = inputUnique[globalPointsNonUnique[i]];
  }

  offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalPoints[i];

  forAll( fluidPatchPointsDispl, i )
  {
    fluidPatchPointsDispl[i] = inputNonUnique[i + offset];
  }

  moveMesh( fluidPatchPointsDispl - fluidPatchPointsDisplOld );
}

void foamFluidSolver::solve(
  const matrix & input,
  matrix & output
  )
{
  fluidPatchPointsDisplOld = fluidPatchPointsDispl;

  vectorField inputNonUnique( globalPointsNonUnique.size(), Foam::vector::zero );

  forAll( inputNonUnique, i )
  {
    for ( int j = 0; j < input.cols(); j++ )
      inputNonUnique[i][j] = input( globalPointsNonUnique[i], j );
  }

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalPoints[i];

  forAll( fluidPatchPointsDispl, i )
  {
    fluidPatchPointsDispl[i] = inputNonUnique[i + offset];
  }

  moveMesh( fluidPatchPointsDispl - fluidPatchPointsDisplOld );

  solve();

  matrix traction;
  getTractionLocal( traction );

  assert( traction.rows() == nGlobalCenters[Pstream::myProcNo()] );

  vectorField writePositionsField( sum( nGlobalCenters ), Foam::vector::zero );

  assert( writePositionsField.size() == data.rows() );

  offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalCenters[i];

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      writePositionsField[i + offset][j] = traction( i, j );

  reduce( writePositionsField, sumOp<vectorField>() );

  output.resize( writePositionsField.size(), 3 );

  for ( int i = 0; i < output.rows(); i++ )
    for ( int j = 0; j < output.cols(); j++ )
      output( i, j ) = writePositionsField[i][j];

  assert( output.rows() == data.rows() );
  assert( output.cols() == data.cols() );

  data = output;
}
