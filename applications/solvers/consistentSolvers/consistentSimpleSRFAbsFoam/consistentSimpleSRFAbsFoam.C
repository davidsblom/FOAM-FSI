/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleSRFFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with single rotating frame.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "SRFModelAbs.H"

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


    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor
            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              //+ SRF->Fcoriolis()
            );
            SRF->addCoriolis(UEqn());

            UEqn().relax();

            solve(UEqn() == -fvc::grad(p));

            //Pressure equation
            p.boundaryField().updateCoeffs();
            volVectorField HU = UEqn().H();
            volScalarField AU = UEqn().A();
            U = HU/AU;
            phi = (fvc::interpolate(HU)/fvc::interpolate(AU)) & mesh.Sf();

            forAll(phi.boundaryField(), patchI)
            {
                if (!phi.boundaryField()[patchI].coupled())
                {
                    phi.boundaryField()[patchI] =
                    (
                        U.boundaryField()[patchI]
                        & mesh.Sf().boundaryField()[patchI]
                    );
                }
            }
            UEqn.clear();

            //Make fluxes relative to SRF frame (and exclude patches which are indicated in dict)
            SRF->relativeFlux(phi);

            adjustPhi(phi, U, p);

            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/fvc::interpolate(AU), p, "laplacian((1|A(U)),p)") == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U -= fvc::grad(p)/AU;
            U.correctBoundaryConditions();
        }


        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
