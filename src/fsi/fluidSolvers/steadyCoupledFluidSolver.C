
/*
 * Class
 *   CoupledFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "steadyCoupledFluidSolver.H"

steadyCoupledFluidSolver::steadyCoupledFluidSolver(
  string name,
  shared_ptr<argList> args,
  shared_ptr<Time> runTime
  )
  :
  CoupledFluidSolver( name, args, runTime )
{}

steadyCoupledFluidSolver::~steadyCoupledFluidSolver()
{}

void steadyCoupledFluidSolver::solve()
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
      fvm::div( phi, U )
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
      underrelaxation = false;
  }

  Info << "Solving for Up, Initial residual = " << initialResidual;
  Info << ", Final residual = " << relativeResidual;
  Info << ", No outer iterations = " << oCorr + 1 << endl;

  continuityErrs();
}
