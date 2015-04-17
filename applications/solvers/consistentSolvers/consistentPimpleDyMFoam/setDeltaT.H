
if ( adjustTimeStep )
{
    scalar maxDeltaTFact = maxCo / (CoNum + SMALL);
    scalar deltaTFact = min( min( maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact ), 1.2 );

    runTime.setDeltaT
    (
        min
        (
            deltaTFact * runTime.deltaT().value(),
            maxDeltaT
        )
    );

    Info << "deltaT = " << runTime.deltaT().value() << endl;
}

// ************************************************************************* //
