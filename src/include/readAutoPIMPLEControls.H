dictionary pimple = mesh.solutionDict().subDict("PIMPLE");

int nOuterCorr(readInt(pimple.lookup("nOuterCorrectors")));//Max outer correctors
int nCorr(readInt(pimple.lookup("nCorrectors")));//Max inner correctors
int nNonOrthCorr = pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);//# non-Orthogonal correctors

scalar momentumTolerance = readScalar(pimple.lookup("momentumTolerance"));//Read in momentumTolerance
scalar pressureTolerance = momentumTolerance;//Default is pressureTolerance equal to momentumTolerance
pimple.readIfPresent("pressureTolerance", pressureTolerance);//Read if present

scalar innerTolerance = 0;//Inner tolerance for inner correctors
pimple.readIfPresent("innerTolerance", innerTolerance);

bool outerLoopConverged=false;//create converged variable for outerloop
scalar innerResidual = 1;//residual initialization for inner correctors
scalar residualPressure = 1;//residual initialization for pressure outer correctors
