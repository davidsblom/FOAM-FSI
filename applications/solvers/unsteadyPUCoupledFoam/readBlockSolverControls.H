
label pRefCell = 0;
scalar pRefValue = 0;
setRefCell( p, mesh.solutionDict().subDict( "blockSolver" ), pRefCell, pRefValue );

scalar convergenceTolerance = readScalar( mesh.solutionDict().subDict( "blockSolver" ).lookup( "convergenceTolerance" ) );
label nOuterCorr = readLabel( mesh.solutionDict().subDict( "blockSolver" ).lookup( "nOuterCorrectors" ) );

// Ensure that the absolute tolerance of the linear solver is less than the
// used convergence tolerance for the non-linear system.
scalar absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "Up" ).lookup( "tolerance" ) );
assert( absTolerance < convergenceTolerance );
