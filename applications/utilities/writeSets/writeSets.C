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
    writePatchPressure

Description
	Write pressure at specified patch with facecentre locations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellSet.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

void readField(word fieldType,word fieldName,const fvMesh& mesh);

int main(int argc, char *argv[])
{	
	timeSelector::addOptions();
#	include "setRootCase.H"	

	// === Create normal mesh and save cellCentres === //
#   include "createTime.H"
	instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

	// === set up directories === //
	fileName setDir("sets");
	if(!isDir(setDir))
	{
		//rmDir(setDir);
		mkDir(setDir);
	}

	// === Read in dictionary === //
    IOdictionary writeSetsDict
    (
        IOobject
        (
            "writeSetsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	wordList fields = writeSetsDict.lookup("fields");
	wordList cellSets = writeSetsDict.lookup("cellSets");
	wordList faceSets = writeSetsDict.lookup("faceSets");
	wordList patchNames = writeSetsDict.lookup("patches");
	
	// === Check if fields are present and from which type === //
	wordList fieldTypes(fields.size());
	wordList fieldNames(fields.size());
	int counter=0;
	forAll(fields,iField){
	    IOobject io
	    (
	        fields[iField],
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::NO_WRITE,
	        true
	    );
	    
	    //Check if is present
	    if (io.headerOk())
        {
            fieldTypes[counter] = io.headerClassName();
            fieldNames[counter] = fields[iField];
            readField(fieldTypes[counter],fieldNames[counter],mesh);
            counter++;
        }
	}
	fieldTypes.setSize(counter);
	fieldNames.setSize(counter);

	
	// === Read in cell sets === //
	wordList cellSetNames(cellSets.size());
	PtrList<cellSet> listCellSet(cellSets.size());
	counter=0;
	forAll(cellSets,iCellSet)
	{
	    IOobject io
	    (
	        cellSets[iCellSet],
	        mesh.time().findInstance
	        (
	        	mesh.dbDir()/mesh.meshSubDir/"sets",
	        	word::null,
	        	IOobject::NO_READ
	       	),
	        mesh.meshSubDir/"sets",
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::NO_WRITE
	    );
		
		if(io.headerOk())
		{
			cellSetNames[counter] == cellSets[iCellSet];
			listCellSet.set(counter, new cellSet(mesh,cellSets[iCellSet],IOobject::MUST_READ));
			counter++;
		}else
		{
			Info << "CellSet " << cellSets[iCellSet] << " not found. " << endl;
		}
	}
	cellSetNames.setSize(counter);
	listCellSet.setSize(counter);
	
	// === Read in face sets === //
	wordList faceSetNames;
	PtrList<faceSet> listFaceSet(faceSets.size());
	counter=0;
	forAll(faceSets,iFaceSet)
	{
	    IOobject io
	    (
	        faceSets[iFaceSet],
	        mesh.time().findInstance
	        (
	        	mesh.dbDir()/mesh.meshSubDir/"sets",
	        	word::null,
	        	IOobject::NO_READ
	       	),
	        mesh.meshSubDir/"sets",
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::NO_WRITE
	    );

		if(io.headerOk())
		{	
			faceSetNames[counter] == faceSets[iFaceSet];
			listFaceSet.set(counter, new faceSet(mesh,faceSets[iFaceSet],IOobject::MUST_READ));
			counter++;
		}else
		{
			Info << "FaceSet " << faceSets[iFaceSet] << " not found. " << endl;
		}
	}
	faceSetNames.setSize(counter);
	listFaceSet.setSize(counter);
	
	// === Get patch IDs === //
	List<label> patchIDs(patchNames.size(),0);
	counter=0;
	forAll(patchNames,iPatch)
	{
		label patchID = mesh.boundaryMesh().findPatchID(patchNames[iPatch]);
		if(patchID>=0){
			patchIDs[counter]=patchID;
			counter++;
		}else
		{
			Info << "Patch " << patchNames[iPatch] << " not found " << endl;
		}
	}
	patchIDs.setSize(counter);	

	// ========================= //
	// === Go over all times === //
	// ========================= //
	forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();
        Info<< "\nTime = " << runTime.timeName() << endl;

		fileName timeDir(setDir+"/"+runTime.timeName());
		if(!isDir(timeDir))
		{
			mkDir(timeDir);
		}
		Info << "time dir created" << endl;

		forAll(fieldNames,iField)
		{
			const word fieldName=fieldNames[iField];
			const word fieldType=fieldTypes[iField];
			if(fieldType=="volScalarField")
			{
				volScalarField field1
				(
					IOobject
					(
						fieldName,
						runTime.timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::NO_WRITE
					),
					mesh
				);
				
				//Write out all cellSets for this field
				for(int i=0;i<listCellSet.size();i++)
				{
					const cellSet& cSet = listCellSet[i];
					if(cSet.size()>0)
					{
						const fileName setFileName(timeDir+"/cellSet_"+cSet.name()+"_"+field1.name()+".xy");
						Info << "Writing " << field1.name() << " for cellSet " << cSet.name() << " in " << setFileName << endl;						
						//If file is already there->delete first
						if(isFile(setFileName)){
							rm(setFileName);
						}
					
						OFstream of(setFileName);
						of << field1.name()<<"\txc\tyc\tzc" << endl;
						forAll(cSet,iCell)
						{
							vector cellCentre = mesh.cellCentres()[cSet.toc()[iCell]];
							scalar FieldCell = field1[cSet.toc()[iCell]];
							of << FieldCell<<"\t"<< cellCentre.x() << "\t" << cellCentre.y() << "\t" << cellCentre.z() << endl;
						}
					}
				}
				
				surfaceScalarField field1Face = fvc::interpolate(field1);
				
				//Write out all faceSets for this field
				for(int i=0;i<listFaceSet.size();i++)
				{
					const faceSet& fSet = listFaceSet[i];
					if(fSet.size()>0)
					{
						const fileName setFileName(timeDir+"/faceSet_"+fSet.name()+"_"+field1.name()+".xy");
						Info << "Writing " << field1.name() << " for faceSet " << fSet.name() << " in " << setFileName << endl;						
						//If file is already there->delete first
						if(isFile(setFileName)){
							rm(setFileName);
						}
						
						OFstream of(setFileName);
						of<<field1.name()<<"\txc\tyc\tzc" << endl;
						forAll(fSet,iFace)
						{
							const label faceID = fSet.toc()[iFace];
							vector faceCentre = mesh.faceCentres()[faceID];
							scalar FieldFace = field1Face[faceID];
							of << FieldFace<<"\t"<< faceCentre.x() << "\t" << faceCentre.y() << "\t" << faceCentre.z() << endl;
						}
					}
				}
				
				forAll(patchIDs,iPatch)
				{
					const label patchID = patchIDs[iPatch];
					const scalarField& fieldPatch(field1.boundaryField()[patchID]);
					const word& patchName = mesh.boundaryMesh()[patchID].name();
					const fileName setFileName(timeDir+"/patch_"+patchName+"_"+field1.name()+".xy");
					Info << "Writing " << field1.name() << " for patch " << patchName << " in " << setFileName << endl;					
					//If file is already there->delete first
					if(isFile(setFileName)){
						rm(setFileName);
					}
					
					//Get facecentres and normals
					const vectorField& faceCentres(mesh.boundaryMesh()[patchID].faceCentres());
					const vectorField& faceNormals(mesh.boundaryMesh()[patchID].faceAreas());
					
					OFstream of(setFileName);
					of<<field1.name()<<"\txc\tyc\tzc\tnx\tny\tnz" << endl;
					forAll(faceCentres,iFace)
					{
						scalar fieldValue=fieldPatch[iFace];
						vector faceNormal=faceNormals[iFace]/mag(faceNormals[iFace]);
						vector faceCentre=faceCentres[iFace];
						of << fieldValue <<"\t" << faceCentre.x() << "\t" << faceCentre.y() << "\t" << faceCentre.z()
						<< "\t"<< faceNormal.x() << "\t"<< faceNormal.y() << "\t"<< faceNormal.z() << endl;
					}
				}
				
			}else if(fieldType=="volVectorField")
			{
				volVectorField field1
				(
					IOobject
					(
						fieldName,
						runTime.timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::NO_WRITE
					),
					mesh
				);
				
				//Write out all cellSets for this field
				for(int i=0;i<listCellSet.size();i++)
				{
					const cellSet& cSet = listCellSet[i];
					if(cSet.size()>0)
					{
						const fileName setFileName(timeDir+"/cellSet_"+cSet.name()+"_"+field1.name()+".xy");
						Info << "Writing " << field1.name() << " for cellSet " << cSet.name() << " in " << setFileName << endl;
						//If file is already there->delete first
						if(isFile(setFileName)){
							rm(setFileName);
						}
					
						OFstream of(setFileName);
						const word& fn = field1.name();
						of << fn<<"x\t"<<fn<<"y\t"<<fn<<"z\txc\tyc\tzc" << endl;
						forAll(cSet,iCell)
						{
							vector cellCentre = mesh.cellCentres()[cSet.toc()[iCell]];
							vector FieldCell = field1[cSet.toc()[iCell]];
							of << FieldCell.x() << "\t" << FieldCell.y() << "\t" << FieldCell.z() <<"\t"<< cellCentre.x() << "\t" << cellCentre.y() << "\t" << cellCentre.z() << endl;
						}
					}
				}
				
				surfaceVectorField field1Face = fvc::interpolate(field1);
				
				//Write out all faceSets for this field
				for(int i=0;i<listFaceSet.size();i++)
				{
					const faceSet& fSet(listFaceSet[i]);
					if(fSet.size()>0)
					{
						const fileName setFileName(timeDir+"/faceSet_"+fSet.name()+"_"+field1.name()+".xy");
						Info << "Writing " << field1.name() << " for faceSet " << fSet.name() << " in " << setFileName << endl;
						//If file is already there->delete first
						if(isFile(setFileName)){
							rm(setFileName);
						}
						
						OFstream of(setFileName);
						const word& fn = field1.name();
						of << fn<<"x\t"<<fn<<"y\t"<<fn<<"z\txc\tyc\tzc" << endl;
						forAll(fSet,iFace)
						{
							const label faceID = fSet.toc()[iFace];
							vector faceCentre = mesh.faceCentres()[faceID];
							vector FieldFace = field1Face[faceID];
							//Info << "faceID(iFace) = " << faceID<<"("<<iFace<<")" << ", faceCentre.x() = " << faceCentre.x() << ", fieldValue = " << FieldFace << endl;
							of << FieldFace.x() << "\t" << FieldFace.y() << "\t" << FieldFace.z() <<"\t"<< faceCentre.x() << "\t" << faceCentre.y() << "\t" << faceCentre.z() << endl;
						}
					}
				}
				
				forAll(patchIDs,iPatch)
				{
					const label patchID = patchIDs[iPatch];
					const vectorField& fieldPatch(field1.boundaryField()[patchID]);
					const word& patchName = mesh.boundaryMesh()[patchID].name();
					const fileName setFileName(timeDir+"/patch_"+patchName+"_"+field1.name()+".xy");
					Info << "Writing " << field1.name() << " for patch " << patchName << " in " << setFileName << endl;
					//If file is already there->delete first
					if(isFile(setFileName)){
						rm(setFileName);
					}
					
					//Get facecentres and normals
					const vectorField& faceCentres(mesh.boundaryMesh()[patchID].faceCentres());
					const vectorField& faceNormals(mesh.boundaryMesh()[patchID].faceAreas());
					
					OFstream of(setFileName);
					const word& fn = field1.name();
					of << fn<<"x\t"<<fn<<"y\t"<<fn<<"z\txc\tyc\tzc\tnx\tny\tnz" << endl;
					forAll(faceCentres,iFace)
					{
						vector fv=fieldPatch[iFace];
						vector faceNormal=faceNormals[iFace]/mag(faceNormals[iFace]);
						vector faceCentre=faceCentres[iFace];
						of << fv.x() << "\t" << fv.y() << "\t" << fv.z()<<"\t" << faceCentre.x() << "\t" << faceCentre.y() << "\t" << faceCentre.z()
						<< "\t"<< faceNormal.x() << "\t"<< faceNormal.y() << "\t"<< faceNormal.z() << endl;
					}
				}
			}
		}	
	}
	
    Info<< "\nEnd\n" << endl;

    return 0;
}

void readField(word fieldType,word fieldName,const fvMesh& mesh){
	if(fieldType=="volScalarField"){
		volScalarField field1
		(
			IOobject
			(
				fieldName,
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE,
				true
			),
			mesh
		);
	}else if(fieldType=="volVectorField"){
		volVectorField field1
		(
			IOobject
			(
				fieldName,
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE,
				true
			),
			mesh
		);
	}
}

// ************************************************************************* //
