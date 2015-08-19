/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    plusPostRANS

Description

    calculates y+ and u+ fields for wall-bounded flows computed with
    one of the available low-Re RANS (no wall function!) turbulence
    models. More specifically it

    :: 	calculates and outputs y+ (avg., min., max.) based on the
	velocity gradient at the wall

    ::	calculates and outputs the wall averaged friction velocity

    ::  writes fields of y+ and U+ to the corresponding time directory

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "nearWallDist.H"
#include "wallDist.H"
#include "wallFvPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"

    scalar matchTolerance = 0;
    scalar matchTol = 1;

    bool noWrite = args.optionFound("noWrite");

    if (args.optionReadIfPresent("tolerance", matchTolerance))
    {
        matchTol = matchTolerance;
        Info<< "Using cell distance to y match tolerance fraction of " << matchTol << nl << endl;
    }else{
        matchTol = 0.001;
        Info<< "Using default cell-to-face distance to y match tolerance fraction of " << matchTol << nl << endl;
    }



    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        wallDist y(mesh, true);

        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" <<endl;
            Info<< "Writing wall distance to field " << y.name() << nl << endl;
            y.write();
        }

		#include "createFields.H"

        const fvPatchList& patches = mesh.boundary();

		Info<< "Summary: " << nl << endl;

        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (isA<wallFvPatch>(currPatch))
            {
                yPlusTemp.boundaryField()[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        RASModel->nu().boundaryField()[patchi]
                       *mag(U.boundaryField()[patchi].snGrad())
                    )
                   /RASModel->nu().boundaryField()[patchi];

				const scalarField& YpTemp = yPlusTemp.boundaryField()[patchi];

                uTau.boundaryField()[patchi] =
                    sqrt
                    (
		    			RASModel->nu()
                   		*mag(U.boundaryField()[patchi].snGrad())
                    );

				const scalarField& uTauTemp = uTau.boundaryField()[patchi];

                /*Info<< "  y+ for Patch " << patchi
                    << " named " << currPatch.name() << ":"
                    << " min: " << min(YpTemp) << " max: " << max(YpTemp)
                    << " average: " << average(YpTemp)
		    		<< nl << endl;

                Info<< "  uTau for Patch " << patchi
                    << " named " << currPatch.name() << ":"
                    << " min: " << min(uTauTemp) << " max: " << max(uTauTemp)
                    << " average: " << average(uTauTemp)
		    		<< nl << endl;*/
            }
        }


		const volVectorField& centers = mesh.C();
		const surfaceVectorField& faceCenters = mesh.Cf();

		//go through all the cells in the mesh
		forAll(uTau, cellI){

			//loop over all the patches
			forAll(patches, patchi){
	            		const fvPatch& currPatch = patches[patchi];

				//loop through all the faces on that patch
				label nFaces = mesh.boundaryMesh()[patchi].size();

				//if this patch is a wall...
				if(isA<wallFvPatch>(currPatch)){
					for(int facei = 0; facei<nFaces; facei++){

					//calculated distance from the current cell to a face on a wall patch
					scalar cellFaceDist ;

					cellFaceDist = Foam::sqrt(sqr(centers[cellI].x()-faceCenters.boundaryField()[patchi][facei].x()) + sqr(centers[cellI].y()-faceCenters.boundaryField()[patchi][facei].y())+ sqr(centers[cellI].z()-faceCenters.boundaryField()[patchi][facei].z()));

					//convert the y value for comparison
					scalar yTemp = y[cellI];

					//fraction difference between our y (i.e. closest perpendicular distance to wall patch) and our seach for the closest wall face
					scalar diffDist = abs(cellFaceDist - yTemp)/max(abs(cellFaceDist),SMALL);

					//if the fraction difference is less than or equal to the match tolerance, search no further.
					if( diffDist <= matchTol){ uTau[cellI] = uTau.boundaryField()[patchi][facei];	break;}

					}//end for loop over faces
				}//end if statement checking isA wall
			}//end loop over patches
		}//end loop over uTau cells

		//true y+ for arbitrary geometry
        yPlus = y.y() * uTau / RASModel->nu();
        //Set yPlus at boundary equal to y of first cell * uTau at wall/nu

        forAll(yPlus.boundaryField(),patchi){
        	const fvPatch& currPatch = patches[patchi];
            if (isA<wallFvPatch>(currPatch))
            {
        		yPlus.boundaryField()[patchi] = y.y().boundaryField()[patchi].patchInternalField() * uTau.boundaryField()[patchi] / RASModel->nu().boundaryField()[patchi];
        	}
        }

        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (isA<wallFvPatch>(currPatch))
            {
		        Info<< "  y+ for Patch " << patchi
		            << " named " << currPatch.name() << ":"
		            << " min: " << min(yPlus.boundaryField()[patchi]) << " max: " << max(yPlus.boundaryField()[patchi])
		            << " average: " << average(yPlus.boundaryField()[patchi])
		    		<< endl;

		        Info<< "  uTau for Patch " << patchi
		            << " named " << currPatch.name() << ":"
					<< " min: " << min(uTau.boundaryField()[patchi]) << " max: " << max(uTau.boundaryField()[patchi])
					<< " average: " << average(uTau.boundaryField()[patchi])
					<< nl << endl;
			}
		}

		//dummy variable holding velocity units
		dimensionedScalar UUnits ( "UUnits", dimensionSet(0,1,-1,0,0), 1.0 );

		//true uPlus over arbitrary geometry
		uPlus = U / stabilise(uTau, SMALL*UUnits);//used to fix divide by zero error if uTau is zero

		if(noWrite){
			Info << "  noWrite option selected, nothing written." << nl <<endl;
		}else{
        	Info << "  Writing yPlus and uPlus to corresponding fields." << nl <<endl;
        	yPlus.write();
        	uPlus.write();
			uTau.write();
		}

    }//end loop over time directories

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
