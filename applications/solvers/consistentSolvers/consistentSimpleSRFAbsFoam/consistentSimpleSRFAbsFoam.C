
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
        #include "initConvergenceCheck.H"

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            p.storePrevIter();
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        #include "convergenceCheck.H"
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
