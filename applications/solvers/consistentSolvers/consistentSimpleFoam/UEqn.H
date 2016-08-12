
// Solve the Momentum equation
tmp<fvVectorMatrix> UEqn
(
    fvm::div( phi, U )
    + turbulence->divDevReff()
);

UEqn().relax();

eqnResidual = solve
    (
    UEqn() == -fvc::grad( p )
    ).initialResidual();

maxResidual = max( eqnResidual, maxResidual );
