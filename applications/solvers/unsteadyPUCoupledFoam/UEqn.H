
// Momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt( U )
    + fvm::div( phi, U )
    + turbulence->divDevReff()
);

UEqn.relax();

UpEqn.insertEquation( 0, UEqn );
