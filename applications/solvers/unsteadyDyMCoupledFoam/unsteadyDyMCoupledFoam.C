
#include "RASModel.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "VectorNFieldTypes.H"
#include "blockLduSolvers.H"
#include "blockVectorNMatrices.H"
#include "dynamicFvMesh.H"
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
  #include "createDynamicFvMesh.H"
  #include "createFields.H"
  #include "initContinuityErrs.H"
  #include "readBlockSolverControls.H"

  Info << "\nStarting time loop\n" << endl;

  while ( runTime.loop() )
  {
    #include "CourantNo.H"

    Info << "Time = " << runTime.timeName() << nl << endl;

    mesh.update();

    #include "volContinuity.H"
    #include "meshCourantNo.H"

    // Outer correction loop to solve the non-linear system
    for ( int oCorr = 0; oCorr < nOuterCorr; oCorr++ )
    {
      // Make the fluxes relative to the mesh motion
      fvc::makeRelative( phi, U );

      #include "calculateCouplingMatrices.H"

      // Initialize block matrix
      #include "initializeBlockMatrix.H"

      p.storePrevIter();

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

      solverPerf.print();

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

      #include "movingMeshContinuityErrs.H"

      // Make the fluxes relative to the mesh motion
      fvc::makeRelative( phi, U );

      turbulence->correct();

      // Make the fluxes absolute
      fvc::makeAbsolute( phi, U );

      if ( currentResidual < convergenceTolerance )
        break;
    }

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }

  Info << "End\n" << endl;

  return 0;
}
