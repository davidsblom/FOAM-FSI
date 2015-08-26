
{
    /*volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) + (turbulence->divDevReff( U ) & U) + fvc::grad( p );
    scalarField magResU = mag( residual.internalField() );
    scalar residualMomentum = std::sqrt( gSumSqr( magResU ) / mesh.globalData().nTotalCells() );

    // == Normalize residual == //
    scalar Uav = ::sqrt( gSumSqr( mag( U.internalField() ) ) / mesh.globalData().nTotalCells() );
    residualMomentum /= ( Uav / runTime.deltaT().value() );*/

    scalar residualMomentum = UEqnPerformance.initialResidual();

    if(oCorr==0)
    {
        initMomentumResidual = residualMomentum;
        initPressureResidual = initInnerResidualPressure;
    }

    // == Check convergence & display residual
    if ( (oCorr > 1 && residualMomentum < max(initMomentumResidual*relativeOuterTolerance, momentumTolerance) && initInnerResidualPressure < max(initPressureResidual*relativeOuterTolerance, pressureTolerance)) || oCorr + 1 == nOuterCorr )
    {
        outerLoopConverged = true;
    }

    Info << "Convergence: momentum residual = " << residualMomentum << "(" << max(initMomentumResidual*relativeOuterTolerance, momentumTolerance) << ") | " << "pressure residual = " << initInnerResidualPressure << "(" << max(initPressureResidual*relativeOuterTolerance, pressureTolerance) << ") | " << "outerCorrectors = (" << oCorr + 1 << "|" << nOuterCorr << "): -> " << outerLoopConverged << endl;
}
