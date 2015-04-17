
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

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

    Info << "\nStarting time loop\n" << endl;

    while ( runTime.run() )
    {
        #include "readTimeControls.H"
        #include "readAutoPIMPLEControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

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
            }
            while ( innerResidual > innerTolerance && corr < nCorr );

            #include "continuityErrs.H"

            // Correct turbulence
            turbulence->correct();

            // Check convergence
            #include "checkPIMPLEResidualConvergence.H"

            oCorr++;
        }
        while ( !outerLoopConverged );

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
