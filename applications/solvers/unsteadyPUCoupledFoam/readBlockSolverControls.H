
label pRefCell = 0;
scalar pRefValue = 0;
setRefCell( p, mesh.solutionDict().subDict( "blockSolver" ), pRefCell, pRefValue );

scalar convergenceTolerance = readScalar( mesh.solutionDict().subDict( "blockSolver" ).lookup( "convergenceTolerance" ) );
label nOuterCorr = readLabel( mesh.solutionDict().subDict( "blockSolver" ).lookup( "nOuterCorrectors" ) );
