
#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "SRFModelAbs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(
    int argc,
    char * argv[]
    )
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for ( runTime++; !runTime.end(); runTime++ )
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "readSIMPLEControls.H"

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor
            tmp<fvVectorMatrix> UEqn
            (
                fvm::div( phi, U )
                + turbulence->divDevReff( U )

                // + SRF->Fcoriolis()
            );
            SRF->addCoriolis( UEqn() );

            UEqn().relax();

            solve( UEqn() == -fvc::grad( p ) );

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
                pEqn.solve();

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
        }


        turbulence->correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
