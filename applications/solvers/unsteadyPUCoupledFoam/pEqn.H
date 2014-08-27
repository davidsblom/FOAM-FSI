
// Pressure parts of the continuity equation
surfaceScalarField rUAf
(
  "rUAf",
  fvc::interpolate( 1.0 / UEqn.A() )
);

surfaceScalarField presSource
(
  "presSource",
  rUAf * ( fvc::interpolate( fvc::grad( p ) ) & mesh.Sf() )
);

fvScalarMatrix pEqn
(
  -fvm::laplacian( rUAf, p )
  ==
  -fvc::div( presSource )
);

pEqn.setReference( pRefCell, pRefValue );

blockMatrixTools::insertEquation( 3, pEqn, A, x, b );
