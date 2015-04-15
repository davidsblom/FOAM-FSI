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
    icoFoam

Description
    Consistent transient solver for incompressible, laminar flow of Newtonian fluids.
    With automatic outer and inner loop corrections.

    By: Thijs Gillebaart, thijsgillebaart@gmail.com


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#   	include "readTimeControls.H"
#   	include "readAutoPISOControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        label oCorr=0;
        do
        {
            Info << "outer iteration: " << oCorr+1 << endl;

            //Solve momentum equation
#           include "UEqn.H"

            // --- PISO loop
            label corr=0;
            do
            {
#               include "pEqn.H"
                corr++;
            }while(innerResidual > innerTolerance && corr < nCorr);for (int corr=0; corr<nCorr; corr++)

            //Check convergence
            #include "checkPISOResidualConvergence.H"

            oCorr++;
        }while(!outerLoopConverged);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
