
tmp<fvVectorMatrix> UEqn
(
    fvm::ddt( U )
    + fvm::div( phi, U )
    + turbulence->divDevReff( U )
);

{ // To ensure S0 and B0 are thrown out of memory
 // Source and boundaryCoeffs need to be saved when relexation is applied to still obtain time consistent behavior.
 // Only source is affected by relaxation, boundaryCoeffs is not relaxatio dependent.
 // BoundaryCoeffs needs to be saved to generate the correct UEqn after solving. Explicit terms (depending on U(n)) need to remain depending on U(n) and not on new solution)
    vectorField S0 = UEqn().source();
    FieldField<Field, vector> B0 = UEqn().boundaryCoeffs();

    UEqn().relax();

    UEqnPerformance = solve( UEqn() == -fvc::grad( p ) );

    // Reset equation to ensure relaxation parameter is not causing problems for time order
    UEqn() =
        (
        fvm::ddt( U )
        + fvm::div( phi, U )
        + turbulence->divDevReff( U )
        );
    UEqn().source() = S0;
    UEqn().boundaryCoeffs() = B0;
}
