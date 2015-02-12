
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
  globalPointsNonUnique( 1, -1 ),
  totalRunTime( 0 ),
  totalNbIterations( 0 )
{
  assert( solidPatchID >= 0 );

  N = getInterfacePointsSize();
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
  displacement.resize( U.boundaryField()[solidPatchID].size(), mesh.nGeometricD() );

  for ( int i = 0; i < displacement.rows(); i++ )
    for ( int j = 0; j < displacement.cols(); j++ )
      displacement( i, j ) = U.boundaryField()[solidPatchID][i][j];
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
  return getInterfaceCentersSize();
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

  readPositions.resize( readPositionsField.size(), mesh.nGeometricD() );

  for ( int i = 0; i < readPositions.rows(); i++ )
    for ( int j = 0; j < readPositions.cols(); j++ )
      readPositions( i, j ) = readPositionsField[i][j];
}

void foamSolidSolver::getReadPositionsLocal( matrix & readPositions )
{
  const vectorField faceCentres( mesh.boundaryMesh()[solidPatchID].faceCentres() );

  readPositions.resize( faceCentres.size(), mesh.nGeometricD() );

  for ( int i = 0; i < readPositions.rows(); i++ )
    for ( int j = 0; j < readPositions.cols(); j++ )
      readPositions( i, j ) = faceCentres[i][j];
}

void foamSolidSolver::getWritePositions( matrix & writePositions )
{
  getReadPositions( writePositions );
}

void foamSolidSolver::getWritePositionsLocal( matrix & writePositions )
{
  getReadPositionsLocal( writePositions );
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
    << " for fluid-solid interface patch, instead "
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
  assert( traction.cols() == mesh.nGeometricD() );

  vectorField tractionSolid( traction.rows(), Foam::vector::zero );

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      tractionSolid[i][j] = traction( i, j );

  tForce->traction() = tractionSolid;
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

  vectorField outputField( getInterfaceCentersSize(), Foam::vector::zero );

  forAll( U.boundaryField()[solidPatchID], i )
  {
    outputField[i + offset] = U.boundaryField()[solidPatchID][i];
  }

  reduce( outputField, sumOp<vectorField>() );

  output.resize( outputField.size(), mesh.nGeometricD() );

  for ( int i = 0; i < output.rows(); i++ )
    for ( int j = 0; j < output.cols(); j++ )
      output( i, j ) = outputField[i][j];

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
