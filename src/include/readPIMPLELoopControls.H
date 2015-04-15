scalar eqnResidual = 1, maxResidual = 0;
scalar outerConvergence = 0;
pimple.readIfPresent("outerConvergence", outerConvergence);
scalar innerConvergence = 0;
pimple.readIfPresent("innerConvergence", innerConvergence);
