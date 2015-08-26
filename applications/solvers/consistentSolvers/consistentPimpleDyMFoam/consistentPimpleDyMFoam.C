
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(
    int argc,
    char * argv[]
    )
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while ( runTime.run() )
    {
        #include "readControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // Make the fluxes absolute
        fvc::makeAbsolute( phi, U );

        // updates V000 if bdf3 is used
        #include "updateV000.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();

        if ( correctPhi && (mesh.moving() || meshChanged) )
        {
            #include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative( phi, U );

        if ( mesh.moving() && checkMeshCourantNo )
        {
            #include "meshCourantNo.H"
        }

        label oCorr = 0;

        do
        {
            Info << "outer iteration: " << oCorr + 1 << endl;

            U.storePrevIter();
            #include "UEqn.H"

            // --- PISO loop
            label corr = 0;

            do
            {
                p.storePrevIter();
                #include "pEqn.H"
                corr++;
            }while(residualPressure > max(initInnerResidualPressure*relativeInnerTolerance, absoluteInnerTolerance) && corr < nCorr);

            #include "movingMeshContinuityErrs.H"

            // Correct turbulence
            turbulence->correct();

            // Check convergence
            #include "checkPIMPLEResidualConvergence.H"

            oCorr++;
        }
        while ( !outerLoopConverged );

        // Update the face velocities
        fvc::makeAbsolute( phi, U );
        #include "updateUf.H"
        fvc::makeRelative( phi, U );

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
