
#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include <assert.h>
#include <yaml-cpp/yaml.h>

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
    #include "readBlockSolverControls.H"

    Info << "\nStarting time loop\n" << endl;

    while ( runTime.loop() )
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        for ( int oCorr = 0; oCorr < nOuterCorr; oCorr++ )
        {
            p.storePrevIter();

            // Initialize the Up block system (matrix, source and reference to Up)
            fvBlockMatrix<vector4> UpEqn( Up );

            // Assemble and insert momentum equation
            #include "UEqn.H"

            // Assemble and insert pressure equation
            #include "pEqn.H"

            // Assemble and insert coupling terms
            #include "couplingTerms.H"

            // Solve the block matrix
            vector4 initialResidual4 = UpEqn.solve().initialResidual();

            // Retrieve solution
            UpEqn.retrieveSolution( 0, U.internalField() );
            UpEqn.retrieveSolution( 3, p.internalField() );

            U.correctBoundaryConditions();
            p.correctBoundaryConditions();

            phi = ( fvc::interpolate( U ) & mesh.Sf() ) + pEqn.flux() + presSource;

            #include "continuityErrs.H"

            p.relax();

            turbulence->correct();

            scalar residual = 0;
            forAll( initialResidual4, i )
            {
                residual = max( residual, initialResidual4[i] );
            }

            if ( residual < convergenceTolerance )
                break;
        }

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}
