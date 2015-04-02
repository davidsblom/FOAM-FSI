
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

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

  Info << "\nStarting time loop\n" << endl;

  while ( runTime.run() )
  {
    #include "readTimeControls.H"
    #include "readPIMPLEControls.H"
    #include "readPIMPLELoopControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"

    runTime++;

    Info << "Time = " << runTime.timeName() << nl << endl;

    label oCorr = 1;
    bool converged = false;
    bool lastIter = false;

    // --- Pressure-velocity PIMPLE corrector loop
    do
    {
      Info << "outer iteration: " << oCorr << endl;

      // reset maxResidual for new outerloop
      maxResidual = 0;

      p.storePrevIter();
      U.storePrevIter();
      #include "UEqn.H"

      // --- PISO loop
      label corr = 0;

      do
      {
        corr++;
        #include "pEqn.H"
      }
      while ( eqnResidual > innerConvergence && corr < nCorr );

      #include "continuityErrs.H"

      turbulence->correct();

      volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) + (turbulence->divDevReff( U ) & U) + fvc::grad( p );

      scalarField magResU = mag( residual.internalField() );
      scalar momentumResidual = std::sqrt( gSumSqr( magResU ) / mesh.globalData().nTotalCells() );
      scalar rmsU = std::sqrt( gSumSqr( mag( U.internalField() ) ) / mesh.globalData().nTotalCells() );
      rmsU /= runTime.deltaT().value();

      // Scale the residual by the root mean square of the velocity field
      momentumResidual /= rmsU;

      // Check if this was the last iter and set converged to true
      if ( lastIter )
      {
        converged = true;
      }

      // Check convergence and set lastIter to true
      if ( momentumResidual < outerConvergence || oCorr == nOuterCorr )
      {
        lastIter = true;
      }

      oCorr++;
    }
    while ( !converged );

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }

  Info << "End\n" << endl;

  return 0;
}
