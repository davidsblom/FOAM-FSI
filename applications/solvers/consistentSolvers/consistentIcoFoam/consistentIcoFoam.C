
#include "fvCFD.H"

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

    while ( runTime.loop() )
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "readAutoPISOControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        label oCorr = 0;

        do
        {
            Info << "outer iteration: " << oCorr + 1 << endl;

            // Solve momentum equation
            #include "UEqn.H"

            // --- PISO loop
            label corr = 0;

            do
            {
                #include "pEqn.H"
                corr++;
            }
            while ( innerResidual > innerTolerance && corr < nCorr );

            for ( int corr = 0; corr < nCorr; corr++ )

            // Check convergence
            #include "checkPISOResidualConvergence.H"

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
