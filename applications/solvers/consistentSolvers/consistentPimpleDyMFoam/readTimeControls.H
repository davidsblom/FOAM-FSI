
Switch adjustTimeStep
(
    runTime.controlDict().lookup( "adjustTimeStep" )
);

scalar maxCo
(
    readScalar( runTime.controlDict().lookup( "maxCo" ) )
);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>( "maxDeltaT", GREAT );

// ************************************************************************* //
