
dictionary pimple = mesh.solutionDict().subDict( "PIMPLE" );

int nOuterCorr( readInt( pimple.lookup( "nOuterCorrectors" ) ) ); // Max outer correctors
int nCorr( readInt( pimple.lookup( "nCorrectors" ) ) ); // Max inner correctors
int nNonOrthCorr = pimple.lookupOrDefault<int>( "nNonOrthogonalCorrectors", 0 ); // # non-Orthogonal correctors

//Outerloop tolerances
scalar momentumTolerance = readScalar( pimple.lookup( "momentumTolerance" ) ); // Read in momentumTolerance
scalar pressureTolerance = momentumTolerance; // Default is pressureTolerance equal to momentumTolerance
pimple.readIfPresent( "pressureTolerance", pressureTolerance ); // Read if present
scalar relativeOuterTolerance = 0.01; // Default is 1e-2
pimple.readIfPresent( "relativeOuterTolerance", relativeOuterTolerance ); // Read if present

//Inner (PISO) loop tolerances
scalar absoluteInnerTolerance = 1e-14; // Inner tolerance for inner correctors
pimple.readIfPresent( "absoluteInnerTolerance", absoluteInnerTolerance );
scalar relativeInnerTolerance = 0.01; // Inner tolerance for inner correctors
pimple.readIfPresent( "relativeInnerTolerance", relativeInnerTolerance );

//initialization
bool outerLoopConverged = false; // create converged variable for outerloop
scalar initInnerResidualPressure = 1;//first pressure residual of outer corrector
scalar residualPressure = 1; // residual initialization for pressure outer correctors
scalar initMomentumResidual = 1;// residual initialization for momentum
scalar initPressureResidual = 1;// residual initialization for pressure outer correctors
label oCorr=0;

lduMatrix::solverPerformance UEqnPerformance;
lduMatrix::solverPerformance pEqnPerformance;
