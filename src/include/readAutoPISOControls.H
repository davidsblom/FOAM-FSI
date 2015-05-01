
dictionary piso = mesh.solutionDict().subDict( "PISO" );

int nOuterCorr( readInt( piso.lookup( "nOuterCorrectors" ) ) ); // Max outer correctors
int nCorr( readInt( piso.lookup( "nCorrectors" ) ) ); // Max inner correctors
int nNonOrthCorr = piso.lookupOrDefault<int>( "nNonOrthogonalCorrectors", 0 ); // # non-Orthogonal correctors

scalar momentumTolerance = readScalar( piso.lookup( "momentumTolerance" ) ); // Read in momentumTolerance
scalar pressureTolerance = momentumTolerance; // Default is pressureTolerance equal to momentumTolerance
piso.readIfPresent( "pressureTolerance", pressureTolerance ); // Read if present

scalar innerTolerance = 0; // Inner tolerance for inner correctors
piso.readIfPresent( "innerTolerance", innerTolerance );

bool outerLoopConverged = false; // create converged variable for outerloop
scalar innerResidual = 1; // residual initialization for inner correctors
scalar residualPressure = 1; // residual initialization for pressure outer correctors
