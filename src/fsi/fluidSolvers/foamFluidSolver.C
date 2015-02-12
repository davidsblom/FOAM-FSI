
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
  fluidPatchPointsDispl( mesh.boundaryMesh()[fluidPatchID].faceCentres().size(), Foam::vector::zero ),
  fluidPatchPointsDisplOld( mesh.boundaryMesh()[fluidPatchID].faceCentres().size(), Foam::vector::zero ),
  nPoints( Pstream::nProcs(), 0 ),
  nGlobalPoints( Pstream::nProcs(), 0 ),
  nGlobalCenters( Pstream::nProcs(), 0 ),
  globalPointsUnique( 1, -1 ),
  globalPointsNonUnique( 1, -1 ),
  totalRunTime( 0 ),
  totalNbIterations( 0 )
{
  assert( fluidPatchID >= 0 );

  N = getInterfaceCentersSize();
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
  return getInterfaceCentersSize();
}

int foamFluidSolver::getInterfacePointsSizeLocal()
{
  return mesh.boundaryMesh()[fluidPatchID].faceCentres().size();
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

  writePositions.resize( writePositionsField.size(), mesh.nGeometricD() );

  for ( int i = 0; i < writePositions.rows(); i++ )
    for ( int j = 0; j < writePositions.cols(); j++ )
      writePositions( i, j ) = writePositionsField[i][j];
}

void foamFluidSolver::getWritePositionsLocal( matrix & writePositions )
{
  const vectorField faceCentres( mesh.boundaryMesh()[fluidPatchID].faceCentres() );

  writePositions.resize( faceCentres.size(), mesh.nGeometricD() );

  for ( int i = 0; i < writePositions.rows(); i++ )
    for ( int j = 0; j < writePositions.cols(); j++ )
      writePositions( i, j ) = faceCentres[i][j];
}

void foamFluidSolver::moveMesh( vectorField motion )
{
  RBFMeshMotionSolver & motionSolver =
    const_cast<RBFMeshMotionSolver &>
    (
    mesh.lookupObject<RBFMeshMotionSolver>( "dynamicMeshDict" )
    );

  Field<vectorField> patches( mesh.boundaryMesh().size(), vectorField( 0 ) );
  patches[fluidPatchID] = motion;
  motionSolver.setMotion( patches );

  mesh.update();
}

void foamFluidSolver::setDisplacementLocal( const matrix & displacement )
{
  fluidPatchPointsDisplOld = fluidPatchPointsDispl;

  for ( int i = 0; i < displacement.rows(); i++ )
    for ( int j = 0; j < displacement.cols(); j++ )
      fluidPatchPointsDispl[i][j] = displacement( i, j );

  moveMesh( fluidPatchPointsDispl - fluidPatchPointsDisplOld );
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

  fluidPatchPointsDisplOld = fluidPatchPointsDispl;

  int offset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    offset += nGlobalCenters[i];

  forAll( fluidPatchPointsDispl, i )
  {
    for ( int j = 0; j < input.cols(); j++ )
      fluidPatchPointsDispl[i][j] = input( i + offset, j );
  }

  moveMesh( fluidPatchPointsDispl - fluidPatchPointsDisplOld );

  solve();

  matrix traction;
  getTractionLocal( traction );

  assert( traction.rows() == nGlobalCenters[Pstream::myProcNo()] );

  vectorField writePositionsField( sum( nGlobalCenters ), Foam::vector::zero );

  assert( writePositionsField.size() == data.rows() );

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      writePositionsField[i + offset][j] = traction( i, j );

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
