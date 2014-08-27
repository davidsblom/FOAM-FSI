
// Momentum equation
fvVectorMatrix UEqn
(
  fvm::ddt( U )
  + fvm::div( phi, U )
  + turbulence->divDevReff( U )
);

UEqn.relax();

blockMatrixTools::insertEquation( 0, UEqn, A, x, b );
