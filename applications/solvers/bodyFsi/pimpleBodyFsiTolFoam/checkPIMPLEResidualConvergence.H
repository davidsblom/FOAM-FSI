
{
    volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) + (turbulence->divDevReff( U ) & U) + fvc::grad( p );
    scalarField magResU = mag( residual.internalField() );
    scalar residualMomentum = std::sqrt( gSumSqr( magResU ) / mesh.globalData().nTotalCells() );

    // == Normalize residual == //
    scalar Uav = ::sqrt( gSumSqr( mag( U.internalField() ) ) / mesh.globalData().nTotalCells() );
    residualMomentum /= ( Uav / runTime.deltaT().value() );

    if(oCorr==0)
    {
        initMomentumResidual = residualMomentum;
        initPressureResidual = initInnerResidualPressure;
    }

    // == Check convergence & display residual
    if ( (oCorr > 1 && residualMomentum < max(initMomentumResidual*relTol,1e-14) && residualPressure < max(initPressureResidual*relTol,1e-14)) || oCorr + 1 == nOuterCorr )
    {
        outerLoopConverged = true;
    }

    Info << "Convergence: momentum residual = " << residualMomentum << "(" << max(initMomentumResidual*relTol,1e-14) << ") | " << "pressure residual = " << residualPressure << "(" << max(initPressureResidual*relTol,1e-14) << ") | " << "outerCorrectors = (" << oCorr + 1 << "|" << nOuterCorr << "): -> " << outerLoopConverged << endl;
}
