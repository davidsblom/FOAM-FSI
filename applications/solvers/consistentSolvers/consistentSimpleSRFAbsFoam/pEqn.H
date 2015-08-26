// Pressure equation
p.boundaryField().updateCoeffs();
volVectorField HU = UEqn().H();
volScalarField AU = UEqn().A();
U = HU / AU;
phi = ( fvc::interpolate( HU ) / fvc::interpolate( AU ) ) & mesh.Sf();

forAll( phi.boundaryField(), patchI )
{
    if ( !phi.boundaryField()[patchI].coupled() )
    {
        phi.boundaryField()[patchI] =
            (
            U.boundaryField()[patchI]
            & mesh.Sf().boundaryField()[patchI]
            );
    }
}
UEqn.clear();

// Make fluxes relative to SRF frame (and exclude patches which are indicated in dict)
SRF->relativeFlux( phi );

adjustPhi( phi, U, p );

// Non-orthogonal pressure corrector loop
for ( int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian( 1.0 / fvc::interpolate( AU ), p, "laplacian((1|A(U)),p)" ) == fvc::div( phi )
    );

    pEqn.setReference( pRefCell, pRefValue );

    // Retain the residual from the first iteration
    if ( nonOrth == 0 )
    {
        eqnResidual = pEqn.solve().initialResidual();
        maxResidual = max( eqnResidual, maxResidual );
    }
    else
    {
        pEqn.solve();
    }

    if ( nonOrth == nNonOrthCorr )
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U -= fvc::grad( p ) / AU;
U.correctBoundaryConditions();
