/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    pimpleFoamRC

Description
    Large time-step transient solver for incompressible, flow using the PIMPLE
    (merged PISO-SIMPLE) algorithm.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "readPIMPLEControls.H"
#       include "readPIMPLELoopControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        label oCorr=1;
        scalar relativeResidual = 1;
        bool converged=false;
        bool lastIter=false;
        // --- Pressure-velocity PIMPLE corrector loop
        do
        {
        	Info << "outer iteration: " << oCorr << endl;
            maxResidual = 0;//reset maxResidual for new outerloop

            p.storePrevIter();
            U.storePrevIter();
#           include "UEqn.H"

            // --- PISO loop
            label corr=0;
            do
            {
                corr++;
#               include "pEqn.H"
            }while(eqnResidual > innerConvergence && corr < nCorr);

#           include "continuityErrs.H"

            turbulence->correct();

            //Get residuals
            scalar residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / (gSumMag( p.internalField())+ SMALL);
            scalar residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / (gSumMag( U.internalField()) + SMALL);
            relativeResidual = max( residualPressure, residualVelocity );
            Info << "residualPressure = " << residualPressure << ", residualVelocity = " << residualVelocity << endl;

            //Check if this was the last iter and set converged to true
            if(lastIter){
                converged = true;
            }
            //Check convergence and set lastIter to true
            if(relativeResidual < outerConvergence || oCorr == nOuterCorr){
                lastIter = true;
            }

            oCorr++;
        }while(!converged);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
