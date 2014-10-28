
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

  // Outer correction loop to solve the non-linear system
  for ( oCorr = 0; oCorr < nOuterCorr; oCorr++ )
  {
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
      fvm::div( phi, U )
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

    continuityErrs();

    turbulence->correct();

    residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / gSumMag( p.internalField() );
    residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / gSumMag( U.internalField() );
    relativeResidual = max( residualPressure, residualVelocity );

    if ( oCorr == 0 )
      initialResidual = relativeResidual;

    if ( relativeResidual < convergenceTolerance )
      break;
  }

  Info << "Solving for Up, Initial residual = " << initialResidual;
  Info << ", Final residual = " << relativeResidual;
  Info << ", No outer iterations = " << oCorr + 1 << endl;
}
