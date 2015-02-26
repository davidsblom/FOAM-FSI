
/*
 * Class
 *   FluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "FluidSolver.H"

FluidSolver::FluidSolver(
  string name,
  std::shared_ptr<argList> args,
  std::shared_ptr<Time> runTime
  )
  :
  foamFluidSolver( name, args, runTime ),
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
  pimple( mesh.solutionDict().subDict( "PIMPLE" ) ),
  nu( transportProperties.lookup( "nu" ) ),
  rho( transportProperties.lookup( "rho" ) ),
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
  AU
  (
  IOobject
  (
    "AU",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
  ),
  mesh,
  1.0 / runTime->deltaT(),
  zeroGradientFvPatchScalarField::typeName
  ),
  HU
  (
  IOobject
  (
    "HU",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
  ),
  mesh,
  U.dimensions() / runTime->deltaT().dimensions(),
  zeroGradientFvPatchVectorField::typeName
  ),
  laminarTransport( U, phi ),
  SfHeader
  (
  "Sf",
  runTime->timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
  Sf
  (
  SfHeader,
  mesh,
  dimensionedVector( "0", dimArea, Foam::vector::zero )
  ),
  correctPhi( false ),
  convergenceTolerance( 1.0e-5 )
{
  initialize();
}

FluidSolver::~FluidSolver()
{}

void FluidSolver::checkTimeDiscretisationScheme()
{
  const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
  const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );

  if ( ddtSchemes.found( "ddt(U)" ) )
  {
    ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
  }
  else
  {
    ddtScheme = word( ddtSchemes.lookup( "default" ) );
  }

  if ( ddtScheme != word( "EulerRC" ) && ddtScheme != word( "backwardRC" ) )
  {
    FatalErrorIn( "checkTimeDiscretisationScheme: " )
    << " ddt(U) scheme is currently set to " << ddtScheme
    << ". This solver only works with ddt(U) scheme = EulerRC or backwardRC."
    << abort( FatalError );
  }
}

void FluidSolver::checkTotalVolume()
{
  scalar newTotalVolume = gSum( mesh.cellVolumes() );
  scalar totalVolRatio = newTotalVolume / totalVolume;

  Info << "Volume: new = " << newTotalVolume << " old = " << totalVolume
       << " change = " << Foam::mag( newTotalVolume - totalVolume )
       << " ratio = " << totalVolRatio - 1 << endl;

  totalVolume = newTotalVolume;
}

void FluidSolver::continuityErrs()
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

void FluidSolver::correctPhiField()
{
  {
    if ( mesh.changing() )
    {
      forAll( U.boundaryField(), patchi )
      {
        if ( U.boundaryField()[patchi].fixesValue() )
        {
          U.boundaryField()[patchi].initEvaluate();
        }
      }

      forAll( U.boundaryField(), patchi )
      {
        if ( U.boundaryField()[patchi].fixesValue() )
        {
          U.correctBoundaryConditions();
          U.boundaryField()[patchi].evaluate();

          phi.boundaryField()[patchi] =
            U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
      }
    }

    wordList pcorrTypes
    (
      p.boundaryField().size(),
      zeroGradientFvPatchScalarField::typeName
    );

    for ( label i = 0; i < p.boundaryField().size(); i++ )
    {
      if ( p.boundaryField()[i].fixesValue() )
      {
        pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
      }
    }

    volScalarField pcorr
    (
      IOobject
      (
        "pcorr",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar( "pcorr", p.dimensions(), 0.0 ),
      pcorrTypes
    );

    for ( int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
    {
      fvScalarMatrix pcorrEqn
      (
        fvm::laplacian( 1.0 / fvc::interpolate( AU ), pcorr, "laplacian((1|A(U)),p)" ) == fvc::div( phi )
      );

      pcorrEqn.setReference( pRefCell, pRefValue );
      pcorrEqn.solve();

      if ( nonOrth == nNonOrthCorr )
      {
        phi -= pcorrEqn.flux();
      }
    }
  }

  continuityErrs();
}

void FluidSolver::courantNo()
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

void FluidSolver::createFields()
{
  U.oldTime();
  phi.oldTime();

  // Creating field to store old face areas needed for ddtPhiCorr
  createSf();

  // Read pressure properties and create turbulence model
  pRefCell = 0;
  pRefValue = 0.0;
  setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );

  turbulence = autoPtr<incompressible::turbulenceModel>
    (
    incompressible::turbulenceModel::New( U, phi, laminarTransport )
    );
}

void FluidSolver::createSf()
{
  if ( !SfHeader.headerOk() )
  {
    const vectorField & allFaceAreas = mesh.faceAreas();

    Sf.internalField() =
      vectorField::subField( allFaceAreas, mesh.nInternalFaces() );

    const fvPatchList & patches = mesh.boundary();

    forAll( patches, patchI )
    {
      Sf.boundaryField()[patchI] =
        patches[patchI].patchSlice( allFaceAreas );
    }
  }
}

void FluidSolver::getAcousticsDensityLocal( matrix & data )
{
  assert( false );
}

void FluidSolver::getAcousticsVelocityLocal( matrix & data )
{
  assert( false );
}

void FluidSolver::getAcousticsPressureLocal( matrix & data )
{
  assert( false );
}

void FluidSolver::getTractionLocal( matrix & traction )
{
  vectorField tractionField( getInterfaceSizeLocal(), Foam::vector::zero );

  int offset = 0;

  forAll( movingPatchIDs, patchI )
  {
    int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

    vectorField tractionFieldPatchI = -rho.value() * nu.value()
      * U.boundaryField()[movingPatchIDs[patchI]].snGrad()
      + rho.value() * p.boundaryField()[patchI]
      * mesh.boundary()[movingPatchIDs[patchI]].nf();

    forAll( tractionFieldPatchI, i )
    {
      tractionField[i + offset] = tractionFieldPatchI[i];
    }

    offset += size;
  }

  assert( tractionField.size() == nGlobalCenters[Pstream::myProcNo()] );

  traction.resize( tractionField.size(), mesh.nGeometricD() );

  for ( int i = 0; i < traction.rows(); i++ )
    for ( int j = 0; j < traction.cols(); j++ )
      traction( i, j ) = tractionField[i][j];
}

void FluidSolver::getWritePositionsLocalAcoustics( matrix & writePositions )
{
  assert( false );
}

void FluidSolver::initContinuityErrs()
{
  sumLocalContErr = 0;
  globalContErr = 0;
  cumulativeContErr = 0;
}

void FluidSolver::initialize()
{
  assert( !init );

  readPIMPLEControls();

  checkTimeDiscretisationScheme();

  createFields();

  readCouplingProperties();

  initContinuityErrs();

  initTotalVolume();

  readTimeControls();
}

void FluidSolver::initTimeStep()
{
  assert( !init );

  timeIndex++;
  t = timeIndex * runTime->deltaT().value();

  Info << "\nTime = " << t << endl;

  readControls();
  checkTotalVolume();
  courantNo();
  setDeltaT();

  init = true;
}

void FluidSolver::initTotalVolume()
{
  totalVolume = sum( mesh.V() ).value();
}

bool FluidSolver::isRunning()
{
  runTime->write();

  Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
       << "  ClockTime = " << runTime->elapsedClockTime() << " s"
       << endl << endl;

  return runTime->loop();
}

void FluidSolver::postProcessing()
{}

void FluidSolver::readControls()
{
  readPIMPLEControls();

  readTimeControls();

  correctPhi = false;

  if ( pimple.found( "correctPhi" ) )
  {
    correctPhi = Switch( pimple.lookup( "correctPhi" ) );
  }
}

void FluidSolver::readCouplingProperties()
{
  if ( couplingProperties.found( "fluidConvergenceTolerance" ) )
    convergenceTolerance = readScalar( couplingProperties.lookup( "fluidConvergenceTolerance" ) );
  else
    FatalErrorIn( "readCouplingProperties" ) << "fluidConvergenceTolerance is not defined" << abort( FatalError );

  // Ensure that the absolute tolerance of the linear solver is less than the
  // used convergence tolerance for the non-linear system.
  scalar absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "U" ).lookup( "tolerance" ) );
  assert( absTolerance < convergenceTolerance );

  absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "p" ).lookup( "tolerance" ) );
  assert( absTolerance < convergenceTolerance );
}

void FluidSolver::readPIMPLEControls()
{
  nOuterCorr = readInt( pimple.lookup( "nOuterCorrectors" ) );
  nCorr = readInt( pimple.lookup( "nCorrectors" ) );

  nNonOrthCorr = pimple.lookupOrDefault<int>( "nNonOrthogonalCorrectors", 0 );

  momentumPredictor = pimple.lookupOrDefault<Switch>( "momentumPredictor", true );

  transonic = pimple.lookupOrDefault<Switch>( "transonic", false );
}

void FluidSolver::readTimeControls()
{
  adjustTimeStep = Switch( runTime->controlDict().lookup( "adjustTimeStep" ) );

  maxCo = readScalar( runTime->controlDict().lookup( "maxCo" ) );

  maxDeltaT = runTime->controlDict().lookupOrDefault<scalar>( "maxDeltaT", GREAT );
}

void FluidSolver::resetSolution()
{
  U == U.oldTime();
  p == p.oldTime();
  phi == phi.oldTime();
}

void FluidSolver::setDeltaT()
{
  if ( adjustTimeStep )
  {
    assert( false );

    scalar maxDeltaTFact = maxCo / (CoNum + SMALL);
    scalar deltaTFact = min( min( maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact ), 1.2 );

    runTime->setDeltaT
    (
      min
      (
        deltaTFact * runTime->deltaT().value(),
        maxDeltaT
      )
    );

    Info << "deltaT = " << runTime->deltaT().value() << endl;
  }
}

void FluidSolver::updateSf()
{
  Sf.oldTime();

  const vectorField & allFaceAreas = mesh.faceAreas();

  Sf.internalField() =
    vectorField::subField( allFaceAreas, mesh.nInternalFaces() );

  const fvPatchList & patches = mesh.boundary();

  forAll( patches, patchI )
  {
    Sf.boundaryField()[patchI] =
      patches[patchI].patchSlice( allFaceAreas );
  }
}

void FluidSolver::solve()
{
  Info << "Solve fluid domain" << endl;

  int ocorr = 0;
  scalar relativeResidual;
  scalar residualPressure;
  scalar residualVelocity;

  lduMatrix::debug = 0;

  // --- PIMPLE loop
  do
  {
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative( phi, U );

    p.storePrevIter();
    U.storePrevIter();

    fvVectorMatrix UEqn
    (
      fvm::ddt( U )
      + fvm::div( phi, U )
      + turbulence->divDevReff( U )
    );

    if ( ocorr != nOuterCorr - 1 )
      UEqn.relax();

    if ( momentumPredictor )
      Foam::solve( UEqn == -fvc::grad( p ) );

    // --- PISO loop
    for ( int corr = 0; corr < nCorr; corr++ )
    {
      HU = UEqn.H();
      AU = UEqn.A();
      U = HU / AU;

      {
        phi = ( fvc::interpolate( HU ) / fvc::interpolate( AU ) ) & mesh.Sf();

        forAll( phi.boundaryField(), patchI )
        {
          if ( !phi.boundaryField()[patchI].coupled() )
          {
            phi.boundaryField()[patchI] =
              (
              U.boundaryField()[patchI]
              & mesh.Sf().boundaryField()[patchI]
              );
          }
        }

        phi += fvc::ddtPhiCorr( 1.0 / AU, U, phi );
      }

      adjustPhi( phi, U, p );

      for ( int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
      {
        fvScalarMatrix pEqn
        (
          fvm::laplacian( 1.0 / fvc::interpolate( AU ), p, "laplacian((1|A(U)),p)" ) == fvc::div( phi )
        );

        pEqn.setReference( pRefCell, pRefValue );

        pEqn.solve( mesh.solutionDict().solver( p.name() ) );

        if ( nonOrth == nNonOrthCorr )
        {
          phi -= pEqn.flux();
        }
      }

      // Explicitly relax pressure for momentum corrector
      if ( ocorr != nOuterCorr - 1 )
      {
        p.relax();
      }

      U -= (1.0 / AU) * fvc::grad( p );
      U.correctBoundaryConditions();
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative( phi, U );

    turbulence->correct();

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute( phi, U );

    residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / (gSumMag( p.internalField() ) + SMALL);
    residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / (gSumMag( U.internalField() ) + SMALL);

    relativeResidual = max( residualPressure, residualVelocity );
  }
  while
  (
    relativeResidual > convergenceTolerance
    && ++ocorr < nOuterCorr
  );

  Info << "Solving for " << U.name()
       << ", Final residual = " << residualVelocity
       << ", No outer iterations " << ocorr + 1 << endl;

  Info << "Solving for " << p.name()
       << ", Final residual = " << residualPressure
       << ", No outer iterations " << ocorr + 1 << endl;

  continuityErrs();

  lduMatrix::debug = 1;
}
