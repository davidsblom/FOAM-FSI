
volVectorField HU = UEqn().H();
volScalarField AU = UEqn().A();
U = HU / AU;

#include "calcPhi.H"

adjustPhi( phi, U, p );

for ( int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian( 1.0 / fvc::interpolate( AU ), p, "laplacian((1|A(U)),p)" ) == fvc::div( phi )
    );

    pEqn.setReference( pRefCell, pRefValue );
    innerResidual = pEqn.solve().initialResidual();

    if ( nonOrth == 0 && corr == 0 )
    {
        residualPressure = innerResidual;
    }

    if ( nonOrth == nNonOrthCorr )
    {
        phi -= pEqn.flux();
    }
}

p.relax();

U -= (1.0 / AU) * fvc::grad( p );
U.correctBoundaryConditions();
