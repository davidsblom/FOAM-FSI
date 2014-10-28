
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "steadySolidSolver.H"

steadySolidSolver::steadySolidSolver (
  string name,
  std::shared_ptr<argList> args,
  std::shared_ptr<Time> runTime
  )
  :
  SolidSolver( name, args, runTime )
{}

steadySolidSolver::~steadySolidSolver()
{}

void steadySolidSolver::solve()
{
  Info << "Solve solid domain" << endl;

  int iCorr = 0;
  scalar relativeResidual = 1;
  scalar initialResidual = 1;
  lduMatrix::debug = 0;

  calculateEpsilonSigma();

  dimensionedVector gravity( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "gravity" ) );

  do
  {
    U.storePrevIter();

    surfaceTensorField shearGradU =
      ( (I - n * n) & fvc::interpolate( gradU ) );

    fvVectorMatrix UEqn
    (
      -fvm::laplacian( 2 * muf + lambdaf, U, "laplacian(DU,U)" )
      ==
      fvc::div(
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

    UEqn.solve();

    U.relax();

    gradU = fvc::grad( U );

    calculateEpsilonSigma();

    relativeResidual = gSumMag( U.internalField() - U.prevIter().internalField() ) / gSumMag( U.internalField() );

    if ( iCorr == 0 )
      initialResidual = relativeResidual;
  }
  while
  (
    relativeResidual > convergenceTolerance
    &&
    ++iCorr < nCorr
  );

  lduMatrix::debug = 1;

  Info << "Solving for " << U.name();
  Info << ", Initial residual = " << initialResidual;
  Info << ", Final residual = " << relativeResidual;
  Info << ", No outer iterations " << iCorr << endl;
}
