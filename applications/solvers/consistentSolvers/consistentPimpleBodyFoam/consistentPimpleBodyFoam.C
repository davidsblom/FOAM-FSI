
/*---------------------------------------------------------------------------*\
*  =========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\    /   O peration     |
\\  /    A nd           | Copyright held by original author
\\\\/     M anipulation  |
\\ -------------------------------------------------------------------------------
\\ License
\\  This file is part of OpenFOAM.
\\
\\  OpenFOAM is free software; you can redistribute it and/or modify it
\\  under the terms of the GNU General Public License as published by the
\\  Free Software Foundation; either version 2 of the License, or (at your
\\  option) any later version.
\\
\\  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
\\  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
\\  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
\\  for more details.
\\
\\  You should have received a copy of the GNU General Public License
\\  along with OpenFOAM; if not, write to the Free Software Foundation,
\\  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\\
\\ Application
\\  pimpleDyMFoam.C
\\
\\ Description
\\  Transient solver for incompressible, flow of Newtonian fluids
\\  on a moving mesh using the PIMPLE (merged PISO-SIMPLE) algorithm.
\\
\\  Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
\\
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"
#include "bodyCollector.H"
#include "RBFMeshMotionSolver.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(
    int argc,
    char * argv[]
    )
{
    #   include "setRootCase.H"
    #   include "createTime.H"
    #   include "createDynamicFvMesh.H"
    #   include "createFields.H"
    #   include "initContinuityErrs.H"
    #   include "readBodyProperties.H" // Reading body properties

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while ( runTime.run() )
    {
        #       include "readControls.H"
        #       include "CourantNo.H"
        #       include "setDeltaT.H"

        fvc::makeAbsolute( phi, U );

        // updates V000 if bdf3 is used
        #include "updateV000.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // Get motion of bodies and update mesh accordingly
        bodyCol.update();
        ms.setMotion( bodyCol.getPatchMotion() );
        bool meshChanged = mesh.update();

        if ( correctPhi && (mesh.moving() || meshChanged) )
        {
            #           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative( phi, U );

        if ( mesh.moving() && checkMeshCourantNo )
        {
            #           include "meshCourantNo.H"
        }

        label oCorr = 0;

        do
        {
            Info << "outer iteration: " << oCorr + 1 << endl;

            U.storePrevIter();
            #           include "UEqn.H"

            // --- PISO loop
            label corr = 0;

            do
            {
                p.storePrevIter();
                #               include "pEqn.H"
                corr++;
            }
            while ( residualPressure > max( initInnerResidualPressure * relativeInnerTolerance, absoluteInnerTolerance ) && corr < nCorr );

            #           include "movingMeshContinuityErrs.H"

            // Correct turbulence
            turbulence->correct();

            // Check convergence
            #           include "checkPIMPLEResidualConvergence.H"

            oCorr++;
        }
        while ( !outerLoopConverged );

        // Update the face velocities
        fvc::makeAbsolute( phi, U );
        #       include "updateUf.H"
        fvc::makeRelative( phi, U );

        runTime.write();

        // Write bodyMotion properties to file (if requested)
        bodyCol.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
