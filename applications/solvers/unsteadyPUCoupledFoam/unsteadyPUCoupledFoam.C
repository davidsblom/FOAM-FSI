
#include "RASModel.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "VectorNFieldTypes.H"
#include "blockLduSolvers.H"
#include "blockVectorNMatrices.H"
#include "volVectorNFields.H"
#include "blockMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(
  int argc,
  char * argv[]
  )
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "createFields.H"
  #include "initContinuityErrs.H"
  #include "readBlockSolverControls.H"

  Info << "\nStarting time loop\n" << endl;

  // Calculate coupling matrices only once since the mesh doesn't change and
  // implicit div and grad operators are only dependant on Sf. Actually
  // coupling terms (div(U) and grad(p)) in blockMatrix do not change,
  // so they could be inserted only once, resetting other parts of
  // blockMatrix to zero at the end of each time step. VV, 30/April/2014
  #include "calculateCouplingMatrices.H"

  while ( runTime.loop() )
  {
    Info << "Time = " << runTime.timeName() << nl << endl;

    lduMatrix::debug = 0;

    int oCorr;

    for ( oCorr = 0; oCorr < nOuterCorr; oCorr++ )
    {
      p.storePrevIter();

      // Initialize block matrix
      #include "initializeBlockMatrix.H"

      // Assemble and insert momentum equation
      #include "UEqn.H"

      // Assemble and insert pressure equation
      #include "pEqn.H"

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

      if ( oCorr == 0 )
      {
        initialResidual = solverPerf.initialResidual();
        initResidual = 0;

        for ( int i = 0; i < 4; i++ )
          initResidual += initialResidual[i];
      }

      currResidual = solverPerf.initialResidual();
      currentResidual = 0;

      for ( int i = 0; i < 4; i++ )
        currentResidual += currResidual[i];

      // Retrieve solution
      blockMatrixTools::retrieveSolution( 0, U.internalField(), Up );
      blockMatrixTools::retrieveSolution( 3, p.internalField(), Up );

      U.correctBoundaryConditions();
      p.correctBoundaryConditions();

      phi = ( fvc::interpolate( U ) & mesh.Sf() ) + pEqn.flux() + presSource;

      p.relax();

      turbulence->correct();

      if ( currentResidual < convergenceTolerance )
        break;
    }

    #include "continuityErrs.H"

    Info << "Solving for Up, Initial residual = " << initialResidual << ", Final residual = " << currResidual << ", No outer iterations = " << oCorr + 1 << endl;

    lduMatrix::debug = 1;

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }

  Info << "End\n" << endl;

  return 0;
}
