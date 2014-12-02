
/*
 * Class
 *   CoupledFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "CoupledFluidSolver.H"

CoupledFluidSolver::CoupledFluidSolver(
  string name,
  shared_ptr<argList> args,
  shared_ptr<Time> runTime
  )
  :
  foamFluidSolver( name, args, runTime ),
  p
  (
  IOobject
  (
    "p",
    runTime->timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
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
  phi
  (
  IOobject
  (
    "phi",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  linearInterpolate( U ) & mesh.Sf()
  ),
  Up
  (
  IOobject
  (
    "Up",
    runTime->timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector4( "zero", dimless, vector4::zero )
  ),
  laminarTransport( U, phi ),
  sumLocalContErr( 0 ),
  globalContErr( 0 ),
  cumulativeContErr( 0 ),
  pRefCell( 0 ),
  pRefValue( 0 ),
  transportProperties
  (
  IOobject
  (
    "transportProperties",
    runTime->constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
  )
  ),
  nu( transportProperties.lookup( "nu" ) ),
  rho( transportProperties.lookup( "rho" ) ),
  convergenceTolerance( readScalar( mesh.solutionDict().subDict( "blockSolver" ).lookup( "convergenceTolerance" ) ) ),
  nOuterCorr( readLabel( mesh.solutionDict().subDict( "blockSolver" ).lookup( "nOuterCorrectors" ) ) ),
  totalVolume( sum( mesh.V() ).value() ),
  CoNum( 0 ),
  meanCoNum( 0 ),
  velMag( 0 )
{
  assert( fluidPatchID >= 0 );

  N = getInterfaceCentersSize();
  dim = 3;
  data.resize( N, dim );
  data.setZero();

  initialize();

  Info << "nu = " << nu << endl;
}

CoupledFluidSolver::~CoupledFluidSolver(){}

void CoupledFluidSolver::checkTotalVolume()
{
  scalar newTotalVolume = gSum( mesh.cellVolumes() );
  scalar totalVolRatio = newTotalVolume / totalVolume;

  Info << "Volume: new = " << newTotalVolume << " old = " << totalVolume
       << " change = " << Foam::mag( newTotalVolume - totalVolume )
       << " ratio = " << totalVolRatio - 1 << endl;

  totalVolume = newTotalVolume;
}

void CoupledFluidSolver::continuityErrs()
{
  volScalarField contErr = fvc::div( phi );

  sumLocalContErr = runTime->deltaT().value() *
    mag( contErr ) ().weightedAverage( mesh.V() ).value();

  globalContErr = runTime->deltaT().value() *
    contErr.weightedAverage( mesh.V() ).value();

  cumulativeContErr += globalContErr;

  Info << "time step continuity errors : sum local = " << sumLocalContErr
       << ", global = " << globalContErr
       << ", cumulative = " << cumulativeContErr
       << endl;
}

void CoupledFluidSolver::courantNo()
{
  CoNum = 0.0;
  meanCoNum = 0.0;
  velMag = 0.0;

  if ( mesh.nInternalFaces() )
  {
    surfaceScalarField magPhi = mag( phi );

    surfaceScalarField SfUfbyDelta =
      mesh.surfaceInterpolation::deltaCoeffs() * magPhi;

    CoNum = max( SfUfbyDelta / mesh.magSf() )
      .value() * runTime->deltaT().value();

    meanCoNum = ( sum( SfUfbyDelta ) / sum( mesh.magSf() ) )
      .value() * runTime->deltaT().value();

    velMag = max( magPhi / mesh.magSf() ).value();
  }

  Info << "Courant Number mean: " << meanCoNum
       << " max: " << CoNum
       << " velocity magnitude: " << velMag
       << endl;
}

void CoupledFluidSolver::createFields()
{
  turbulence = autoPtr<incompressible::turbulenceModel>
    (
    incompressible::turbulenceModel::New( U, phi, laminarTransport )
    );
}

void CoupledFluidSolver::initialize()
{
  createFields();
}

void CoupledFluidSolver::getTractionLocal( matrix & traction )
{
  vectorField tractionField = -rho.value() * nu.value()
    * U.boundaryField()[fluidPatchID].snGrad()
    + rho.value() * p.boundaryField()[fluidPatchID]
    * mesh.boundary()[fluidPatchID].nf();

  assert( tractionField.size() == nGlobalCenters[Pstream::myProcNo()] );

  traction.resize( tractionField.size(), 3 );

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      traction( i, j ) = tractionField[i][j];
}

void CoupledFluidSolver::initTimeStep()
{
  assert( !init );

  timeIndex++;
  t = timeIndex * runTime->deltaT().value();

  Info << "\nTime = " << t << endl;

  checkTotalVolume();
  courantNo();

  init = true;
}

bool CoupledFluidSolver::isRunning()
{
  runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;
         
  return runTime->loop();
}

void CoupledFluidSolver::moveMesh( vectorField motion )
{
  bool fvMotionSolver =
    mesh.objectRegistry::foundObject<pointVectorField>
    (
    "pointMotionU"
    );

  bool feMotionSolver =
    mesh.objectRegistry::foundObject<tetPointVectorField>
    (
    "motionU"
    );

  bool rbfMotionSolver =
    mesh.objectRegistry::foundObject<RBFMotionSolverExt>( "dynamicMeshDict" );

  bool idwMotionSolver =
    mesh.objectRegistry::foundObject<IDWMotionSolver>( "dynamicMeshDict" );

  if ( fvMotionSolver )
  {
    assert( !rbfMotionSolver );
    assert( !feMotionSolver );
    assert( !idwMotionSolver );

    // Move whole fluid mesh
    pointField newPoints = mesh.allPoints();

    const labelList & meshPoints = mesh.boundaryMesh()[fluidPatchID].meshPoints();

    forAll( motion, pointI )
    {
      newPoints[meshPoints[pointI]] -= motion[pointI];
    }

    twoDPointCorrector twoDCorrector( mesh );

    twoDCorrector.correctPoints( newPoints );

    mesh.movePoints( newPoints );

    pointVectorField & motionU =
      const_cast<pointVectorField &>
      (
      mesh.objectRegistry::
      lookupObject<pointVectorField>
      (
        "pointMotionU"
      )
      );

    fixedValuePointPatchVectorField & motionUFluidPatch =
      refCast<fixedValuePointPatchVectorField>
      (
      motionU.boundaryField()[fluidPatchID]
      );

    motionUFluidPatch == motion / runTime->deltaT().value();
  }

  if ( feMotionSolver )
  {
    assert( !rbfMotionSolver );
    assert( !fvMotionSolver );
    assert( !idwMotionSolver );

    // Move whole fluid mesh
    pointField newPoints = mesh.allPoints();

    const labelList & meshPoints = mesh.boundaryMesh()[fluidPatchID].meshPoints();

    forAll( motion, pointI )
    {
      newPoints[meshPoints[pointI]] -= motion[pointI];
    }

    twoDPointCorrector twoDCorrector( mesh );

    twoDCorrector.correctPoints( newPoints );

    mesh.movePoints( newPoints );

    tetPointVectorField & motionU =
      const_cast<tetPointVectorField &>
      (
      mesh.objectRegistry::
      lookupObject<tetPointVectorField>
      (
        "motionU"
      )
      );

    fixedValueTetPolyPatchVectorField & motionUFluidPatch =
      refCast<fixedValueTetPolyPatchVectorField>
      (
      motionU.boundaryField()[fluidPatchID]
      );

    tetPolyPatchInterpolation tppi
    (
      refCast<const faceTetPolyPatch>( motionUFluidPatch.patch() )
    );

    motionUFluidPatch ==
    tppi.pointToPointInterpolate
    (
      motion / runTime->deltaT().value()
    );
  }

  if ( rbfMotionSolver )
  {
    assert( !fvMotionSolver );
    assert( !feMotionSolver );
    assert( !idwMotionSolver );

    RBFMotionSolverExt & motionSolver =
      const_cast<RBFMotionSolverExt &>
      (
      mesh.lookupObject<RBFMotionSolverExt>( "dynamicMeshDict" )
      );

    Field<vectorField> patches( mesh.boundaryMesh().size(), vectorField( 0 ) );
    patches[fluidPatchID] = motion;
    motionSolver.setMotion( patches );
  }

  if ( idwMotionSolver )
  {
    assert( !fvMotionSolver );
    assert( !feMotionSolver );
    assert( !rbfMotionSolver );

    IDWMotionSolver & motionSolver =
      const_cast<IDWMotionSolver &>
      (
      mesh.lookupObject<IDWMotionSolver>( "dynamicMeshDict" )
      );

    Field<vectorField> patches( mesh.boundaryMesh().size(), vectorField( 0 ) );
    patches[fluidPatchID] = motion;
    motionSolver.setMotion( patches );
  }

  assert( rbfMotionSolver || fvMotionSolver || feMotionSolver || idwMotionSolver );

  mesh.update();
}

void CoupledFluidSolver::readBlockSolverControls()
{
  setRefCell( p, mesh.solutionDict().subDict( "blockSolver" ), pRefCell, pRefValue );
}

void CoupledFluidSolver::resetSolution()
{
  U == U.oldTime();
  p == p.oldTime();
  phi == phi.oldTime();
}

void CoupledFluidSolver::solve()
{
  Info << "Solve fluid domain" << endl;

  int oCorr;
  scalar initialResidual = 1;
  scalar relativeResidual = 1;
  scalar residualPressure = 1;
  scalar residualVelocity = 1;

  // Outer correction loop to solve the non-linear system
  for ( oCorr = 0; oCorr < nOuterCorr; oCorr++ )
  {
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative( phi, U );

    blockVectorMatrix pInU( fvm::grad( p ) );
    blockVectorMatrix UInp( fvm::div( U ) );

    // Initialize block matrix

    // Create block matrix
    BlockLduMatrix<vector4> A( mesh );

    // Block matrix - create x and b
    vector4Field & x = Up.internalField();
    vector4Field b( x.size(), vector4::zero );

    // Set block interfaces properly
    A.interfaces() = Up.boundaryField().blockInterfaces();

    p.storePrevIter();
    U.storePrevIter();

    // Assemble and insert momentum equation

    // Momentum equation
    fvVectorMatrix UEqn
    (
      fvm::ddt( U )
      + fvm::div( phi, U )
      + turbulence->divDevReff( U )
    );

    UEqn.relax();

    blockMatrixTools::insertEquation( 0, UEqn, A, x, b );

    // Assemble and insert pressure equation

    // Pressure parts of the continuity equation
    surfaceScalarField rUAf
    (
      "rUAf",
      fvc::interpolate( 1.0 / UEqn.A() )
    );

    surfaceScalarField presSource
    (
      "presSource",
      rUAf * ( fvc::interpolate( fvc::grad( p ) ) & mesh.Sf() )
    );

    fvScalarMatrix pEqn
    (
      -fvm::laplacian( rUAf, p )
      ==
      -fvc::div( presSource )
    );

    pEqn.setReference( pRefCell, pRefValue );

    blockMatrixTools::insertEquation( 3, pEqn, A, x, b );

    // Insert coupling, updating the boundary contributions
    // Last argument in insertBlockCoupling says if the first location
    // should be incremented. This is needed for arbitrary positioning
    // of U and p in the system. This could be better. VV, 30/April/2014
    blockMatrixTools::insertBlockCoupling( 3, 0, UInp, U, A, b, false );
    blockMatrixTools::insertBlockCoupling( 0, 3, pInU, p, A, b, true );

    // Solve the block matrix
    BlockSolverPerformance<vector4> solverPerf =
      BlockLduSolver<vector4>::New
      (
      word( "Up" ),
      A,
      mesh.solutionDict().solver( "Up" )
      )->solve( Up, b );
    solverPerf.print();

    // Retrieve solution
    blockMatrixTools::retrieveSolution( 0, U.internalField(), Up );
    blockMatrixTools::retrieveSolution( 3, p.internalField(), Up );

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = ( fvc::interpolate( U ) & mesh.Sf() ) + pEqn.flux() + presSource;

    p.relax();

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative( phi, U );

    turbulence->correct();

    // Make the fluxes absolute
    fvc::makeAbsolute( phi, U );

    residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / (gSumMag( p.internalField() ) + SMALL);
    residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / (gSumMag( U.internalField() ) + SMALL);
    relativeResidual = max( residualPressure, residualVelocity );

    scalar initResidual = 0;
    forAll( solverPerf.initialResidual(), i )
    {
      initResidual += solverPerf.initialResidual()[i];
    }

    relativeResidual = min( relativeResidual, initResidual );

    if ( oCorr == 0 )
      initialResidual = relativeResidual;

    if ( relativeResidual < convergenceTolerance && oCorr > 1 )
      break;
  }

  Info << "Solving for Up, Initial residual = " << initialResidual;
  Info << ", Final residual = " << relativeResidual;
  Info << ", No outer iterations = " << oCorr + 1 << endl;

  continuityErrs();
}
