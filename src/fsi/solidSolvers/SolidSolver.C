
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SolidSolver.H"

SolidSolver::SolidSolver (
  string name,
  std::shared_ptr<argList> args,
  std::shared_ptr<Time> runTime
  )
  :
  foamSolidSolver( name, args, runTime ),
  gradU( fvc::grad( U ) ),
  epsilon
  (
  IOobject
  (
    "epsilonGreen",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedSymmTensor( "zero", dimless, symmTensor::zero )
  ),
  sigma
  (
  IOobject
  (
    "sigma",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedSymmTensor( "zero", dimForce / dimArea, symmTensor::zero )
  ),
  rheology( sigma, U ),
  rho( rheology.rho() ),
  mu( rheology.mu() ),
  lambda( rheology.lambda() ),
  muf( fvc::interpolate( mu, "mu" ) ),
  lambdaf( fvc::interpolate( lambda, "lambda" ) ),
  n( mesh.Sf() / mesh.magSf() ),
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
  )
{}

SolidSolver::~SolidSolver()
{}

void SolidSolver::calculateEpsilonSigma()
{
  // - Green finite strain tensor
  epsilon = symm( gradU ) + 0.5 * symm( gradU & gradU.T() );

  // - second Piola-Kirchhoff stress tensor
  sigma = 2 * mu * epsilon + lambda * ( I * tr( epsilon ) );
}

void SolidSolver::initialize()
{
  readCouplingProperties();
}

void SolidSolver::initTimeStep()
{
  assert( !init );

  readSolidMechanicsControls();

  init = true;
}

bool SolidSolver::isRunning()
{
  return runTime->loop();
}

void SolidSolver::readSolidMechanicsControls()
{
  const dictionary & stressControl =
    mesh.solutionDict().subDict( "solidMechanics" );

  nCorr = readInt( stressControl.lookup( "nCorrectors" ) );
  convergenceTolerance = readScalar( stressControl.lookup( "U" ) );
}

void SolidSolver::solve()
{
  Info << "Solve solid domain" << endl;

  int iCorr = 0;
  scalar initialResidual = 0;
  lduMatrix::solverPerformance solverPerf;
  lduMatrix::debug = 0;

  U == U.oldTime();
  gradU = fvc::grad( U );

  calculateEpsilonSigma();

  dimensionedVector gravity( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "gravity" ) );

  do
  {
    U.storePrevIter();

    surfaceTensorField shearGradU =
      ( (I - n * n) & fvc::interpolate( gradU ) );

    fvVectorMatrix UEqn
    (
      rho * fvm::d2dt2( U )
      ==
      fvm::laplacian( 2 * muf + lambdaf, U, "laplacian(DU,U)" )
      + fvc::div(
        mesh.magSf()
        * (
          -(muf + lambdaf) * ( fvc::snGrad( U ) & (I - n * n) )
          + lambdaf * tr( shearGradU & (I - n * n) ) * n
          + muf * (shearGradU & n)
          + muf * ( n & fvc::interpolate( gradU & gradU.T() ) )
          + 0.5 * lambdaf
          * ( n * tr( fvc::interpolate( gradU & gradU.T() ) ) )
          + ( n & fvc::interpolate( sigma & gradU ) )
          )
        )
    );

    // Add gravity

    UEqn -= rho * gravity;

    solverPerf = UEqn.solve();

    if ( iCorr == 0 )
    {
      initialResidual = solverPerf.initialResidual();
    }

    U.relax();

    gradU = fvc::grad( U );

    calculateEpsilonSigma();
  }
  while
  (
    solverPerf.initialResidual() > convergenceTolerance
    &&
    ++iCorr < nCorr
  );

  lduMatrix::debug = 1;

  Info << "Solving for " << U.name()
       << ", Initial residual = " << initialResidual
       << ", Final residual = " << solverPerf.initialResidual()
       << ", No outer iterations " << iCorr << endl;
}
