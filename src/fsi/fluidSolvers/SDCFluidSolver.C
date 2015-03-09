
/*
 * Class
 *   SDCFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCFluidSolver.H"

SDCFluidSolver::SDCFluidSolver(
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
  Uf
  (
  IOobject
  (
    "Uf",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  fvc::interpolate( U )
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
  rhsU
  (
  IOobject
  (
    "rhsU",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector( "rhsU", dimVelocity / dimTime, Foam::vector::zero ),
  zeroGradientFvPatchVectorField::typeName
  ),
  rhsUf
  (
  IOobject
  (
    "rhsUf",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  fvc::interpolate( rhsU )
  ),
  convergenceTolerance( 1.0e-5 ),
  k( 0 ),
  pStages(),
  phiStages(),
  UStages(),
  UfStages()
{
  initialize();
}

SDCFluidSolver::~SDCFluidSolver()
{}

void SDCFluidSolver::createFields()
{
  U.oldTime();
  phi.oldTime();

  // Read pressure properties and create turbulence model
  pRefCell = 0;
  pRefValue = 0.0;
  setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );
}

void SDCFluidSolver::getAcousticsDensityLocal( matrix & data )
{
  assert( false );
}

void SDCFluidSolver::getAcousticsVelocityLocal( matrix & data )
{
  assert( false );
}

void SDCFluidSolver::getAcousticsPressureLocal( matrix & data )
{
  assert( false );
}

void SDCFluidSolver::getTractionLocal( matrix & traction )
{
  assert( false );
}

void SDCFluidSolver::getWritePositionsLocalAcoustics( matrix & writePositions )
{
  assert( false );
}

void SDCFluidSolver::initialize()
{
  assert( !init );

  readPIMPLEControls();

  createFields();

  readCouplingProperties();
}

void SDCFluidSolver::initTimeStep()
{
  assert( !init );

  readPIMPLEControls();

  init = true;
}

bool SDCFluidSolver::isRunning()
{
  runTime->write();

  Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
       << "  ClockTime = " << runTime->elapsedClockTime() << " s"
       << endl << endl;

  return runTime->loop();
}

void SDCFluidSolver::readCouplingProperties()
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

void SDCFluidSolver::readPIMPLEControls()
{
  nOuterCorr = readInt( pimple.lookup( "nOuterCorrectors" ) );
  nCorr = readInt( pimple.lookup( "nCorrectors" ) );

  nNonOrthCorr = pimple.lookupOrDefault<int>( "nNonOrthogonalCorrectors", 0 );
}

void SDCFluidSolver::resetSolution()
{
  assert( false );
}

void SDCFluidSolver::setNumberOfStages( int k )
{
  this->k = k;

  for ( int i = 0; i < k; i++ )
  {
    pStages.push_back( volScalarField( p ) );
    phiStages.push_back( surfaceScalarField( phi ) );
    UStages.push_back( volVectorField( U ) );
    UfStages.push_back( surfaceVectorField( Uf ) );
  }
}

void SDCFluidSolver::nextTimeStep()
{
  timeIndex++;

  if ( pStages.size() == static_cast<unsigned>(k) )
  {
    for ( int i = 0; i < k; i++ )
    {
      pStages.at( i ) = p;
      phiStages.at( i ) = phi;
      UStages.at( i ) = U;
      UfStages.at( i ) = Uf;
    }
  }
}

void SDCFluidSolver::solve()
{
  assert( false );
}

void SDCFluidSolver::finalizeTimeStep()
{
  foamFluidSolver::finalizeTimeStep();
}

int SDCFluidSolver::getDOF()
{
  return U.size() * 3 + Uf.size() * 3;
}

void SDCFluidSolver::getSolution( Eigen::VectorXd & solution )
{
  for ( int i = 0; i < U.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      solution( i * 3 + j ) = U[i][j];

  for ( int i = 0; i < Uf.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      solution( i * 3 + j + U.size() * 3 ) = Uf[i][j];
}

double SDCFluidSolver::getTimeStep()
{
  return runTime->deltaT().value();
}

void SDCFluidSolver::evaluateFunction(
  const int k,
  const Eigen::VectorXd & q,
  const double t,
  Eigen::VectorXd & f
  )
{
  p = pStages.at( k );
  phi = phiStages.at( k );

  for ( int i = 0; i < U.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      U[i][j] = q( i * 3 + j );

  for ( int i = 0; i < Uf.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      Uf[i][j] = q( i * 3 + j + U.size() * 3 );

  volVectorField F = fvc::laplacian( nu, U ) - fvc::div( phi, U ) - fvc::grad( p );
  surfaceVectorField Ff = fvc::interpolate( F );

  for ( int i = 0; i < F.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      f( i * 3 + j ) = F[i][j];

  for ( int i = 0; i < Ff.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      f( i * 3 + j + F.size() * 3 ) = Ff[i][j];
}

void SDCFluidSolver::implicitSolve(
  const int k,
  const double t,
  const double dt,
  const Eigen::VectorXd & qold,
  const Eigen::VectorXd & rhs,
  Eigen::VectorXd & f,
  Eigen::VectorXd & result
  )
{
  runTime->setDeltaT( dt );
  runTime->setTime( t, timeIndex );

  p = pStages.at( k + 1 );
  phi = phiStages.at( k + 1 );
  U = UStages.at( k + 1 );
  Uf = UfStages.at( k + 1 );

  for ( int i = 0; i < U.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      U.oldTime()[i][j] = qold( i * 3 + j );

  for ( int i = 0; i < Uf.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      Uf.oldTime()[i][j] = qold( i * 3 + j + U.size() * 3 );

  for ( int i = 0; i < U.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      rhsU[i][j] = rhs( i * 3 + j );

  for ( int i = 0; i < Uf.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      rhsUf[i][j] = rhs( i * 3 + j + U.size() * 3 );

  // PIMPLE algorithm

  for ( label oCorr = 0; oCorr < nOuterCorr; oCorr++ )
  {
    // for computing Hp
    fvVectorMatrix UEqn
    (
      fvm::div( phi, U )
      - fvm::laplacian( nu, U )
    );

    // for computing \tilde{a}_{P}
    fvVectorMatrix UEqnt
    (
      fvm::ddt( U )
      + fvm::div( phi, U )
      - fvm::laplacian( nu, U )
    );

    // Eqn 25
    Foam::solve( UEqnt == -fvc::grad( p ) + rhsU / dt );

    // --- PISO loop

    for ( label corr = 0; corr < nCorr; corr++ )
    {
      volVectorField HU = UEqn.H();       // Eqn 6
      volScalarField AU = UEqnt.A();       // Eqn 14

      dimensionedScalar rDeltaT = 1.0 / runTime->deltaT();

      // Eqn 15
      U = HU / AU + rDeltaT / AU * U.oldTime() + 1 / AU * rhsU / dt;

      // Eqn 16
      Uf = fvc::interpolate( HU ) / fvc::interpolate( AU ) + rDeltaT / fvc::interpolate( AU ) * Uf.oldTime() + 1 / fvc::interpolate( AU ) * rhsUf / dt;

      {
        forAll( Uf.boundaryField(), patchI )
        {
          if ( !Uf.boundaryField()[patchI].coupled() )
          {
            Uf.boundaryField()[patchI] =
              (
              U.boundaryField()[patchI]
              );
          }
        }

        phi = Uf & mesh.Sf();

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
      }

      for ( label nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
      {
        fvScalarMatrix pEqn
        (
          fvm::laplacian
          (
            1.0 / fvc::interpolate( AU ), p,
            "laplacian((1|A(U)),p)"
          )
          == fvc::div( phi )
        );

        pEqn.setReference( pRefCell, pRefValue );
        pEqn.solve();

        if ( nonOrth == nNonOrthCorr )
        {
          // Eqn 28
          Uf -= pEqn.flux() * mesh.Sf() / ( mesh.magSf() * mesh.magSf() );

          // Eqn 29
          phi = Uf & mesh.Sf();
        }
      }

      // Eqn 30
      U -= fvc::grad( p ) / AU;

      U.correctBoundaryConditions();
    }
  }

  pStages.at( k + 1 ) = p;
  phiStages.at( k + 1 ) = phi;
  UStages.at( k + 1 ) = U;
  UfStages.at( k + 1 ) = Uf;

  for ( int i = 0; i < U.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      result( i * 3 + j ) = U[i][j];

  for ( int i = 0; i < Uf.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      result( i * 3 + j + U.size() * 3 ) = Uf[i][j];

  volVectorField F = fvc::laplacian( nu, U ) - fvc::div( phi, U ) - fvc::grad( p );
  dimensionedScalar rDeltaT = 1.0 / runTime->deltaT();
  surfaceVectorField Ff = rDeltaT * ( Uf - Uf.oldTime() ) - rhsUf / dt;

  for ( int i = 0; i < F.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      f( i * 3 + j ) = F[i][j];

  for ( int i = 0; i < Ff.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      f( i * 3 + j + F.size() * 3 ) = Ff[i][j];
}
