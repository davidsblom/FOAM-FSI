
label pRefCell = 0;
scalar pRefValue = 0;
setRefCell( p, mesh.solutionDict().subDict( "blockSolver" ), pRefCell, pRefValue );

scalar convergenceTolerance = readScalar( mesh.solutionDict().subDict( "blockSolver" ).lookup( "convergenceTolerance" ) );
label nOuterCorr = readLabel( mesh.solutionDict().subDict( "blockSolver" ).lookup( "nOuterCorrectors" ) );

scalar initResidual = 1;
scalar currentResidual = 1;
Foam::VectorN<double, 4> initialResidual;
Foam::VectorN<double, 4> currResidual;
