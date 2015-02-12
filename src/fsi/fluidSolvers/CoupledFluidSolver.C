
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
    IOobject::NO_WRITE
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
    IOobject::NO_WRITE
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

  initialize();
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

void CoupledFluidSolver::getAcousticsDensityLocal( matrix & data )
{
  assert( false );
}

void CoupledFluidSolver::getAcousticsVelocityLocal( matrix & data )
{
  assert( false );
}

void CoupledFluidSolver::getAcousticsPressureLocal( matrix & data )
{
  assert( false );
}

void CoupledFluidSolver::getTractionLocal( matrix & traction )
{
  vectorField tractionField = -rho.value() * nu.value()
    * U.boundaryField()[fluidPatchID].snGrad()
    + rho.value() * p.boundaryField()[fluidPatchID]
    * mesh.boundary()[fluidPatchID].nf();

  assert( tractionField.size() == nGlobalCenters[Pstream::myProcNo()] );

  traction.resize( tractionField.size(), mesh.nGeometricD() );

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      traction( i, j ) = tractionField[i][j];
}

void CoupledFluidSolver::getWritePositionsLocalAcoustics( matrix & writePositions )
{
  assert( false );
}

void CoupledFluidSolver::initTimeStep()
{
  assert( !init );

  timeIndex++;
  t = timeIndex * runTime->deltaT().value();

  Info << "\nTime = " << runTime->timeName() << endl;

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
  bool underrelaxation = true;

  // Outer correction loop to solve the non-linear system
  for ( oCorr = 0; oCorr < nOuterCorr; oCorr++ )
  {
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative( phi, U );

    p.storePrevIter();
    U.storePrevIter();

    // Initialize the Up block system (matrix, source and reference to Up)
    fvBlockMatrix<vector4> UpEqn( Up );

    // Momentum equation
    fvVectorMatrix UEqn
    (
      fvm::ddt( U )
      + fvm::div( phi, U )
      + turbulence->divDevReff( U )
    );

    if ( underrelaxation )
      UEqn.relax();

    UpEqn.insertEquation( 0, UEqn );

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
      rUAf * ( fvc::interpolate( fvc::grad( p, "grad(pSource)" ) ) & mesh.Sf() )
    );

    fvScalarMatrix pEqn
    (
      -fvm::laplacian( rUAf, p )
      ==
      -fvc::div( presSource )
    );

    pEqn.setReference( pRefCell, pRefValue );

    UpEqn.insertEquation( 3, pEqn );

    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
    BlockLduSystem<Foam::vector, Foam::vector> pInU( fvm::grad( p ) );

    // Calculate div U coupling.  Could be calculated only once since
    // it is only geometry dependent.  VV, 9/June/2014
    BlockLduSystem<Foam::vector, scalar> UInp( fvm::UDiv( U ) );

    // Last argument in insertBlockCoupling says if the column direction
    // should be incremented. This is needed for arbitrary positioning
    // of U and p in the system. This could be better. VV, 30/April/2014
    UpEqn.insertBlockCoupling( 0, 3, pInU, true );
    UpEqn.insertBlockCoupling( 3, 0, UInp, false );

    // Solve the block matrix
    vector4 initialResidual4 = UpEqn.solve().initialResidual();

    // UpEqn.initialResidual();

    // Retrieve solution
    UpEqn.retrieveSolution( 0, U.internalField() );
    UpEqn.retrieveSolution( 3, p.internalField() );

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = ( fvc::interpolate( U ) & mesh.Sf() ) + pEqn.flux() + presSource;

    if ( underrelaxation )
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
    forAll( initialResidual4, i )
    {
      initResidual = max( initResidual, initialResidual4[i] );
    }

    relativeResidual = min( relativeResidual, initResidual );

    if ( oCorr == 0 )
      initialResidual = relativeResidual;

    if ( relativeResidual < convergenceTolerance && oCorr > 1 && !underrelaxation )
      break;

    if ( relativeResidual < convergenceTolerance && oCorr > 1 )
      break;
  }

  Info << "Solving for Up, Initial residual = " << initialResidual;
  Info << ", Final residual = " << relativeResidual;
  Info << ", No outer iterations = " << oCorr + 1 << endl;

  continuityErrs();
}
