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
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
bodyCollector

\*---------------------------------------------------------------------------*/

#include "bodyCollector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(bodyCollector, 0);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void bodyCollector::setAllMotion(Field<vectorField>& allMotionPerPatch,const Field<vectorField>& motionPerPatch) const{

    if(motionPerPatch.size() != mesh_.boundaryMesh().size()){
        FatalErrorIn("tmp<Field<vectorField> > bodyCollector::setAllMotion(Field<vectorField>& allMotionPerPatch,const Field<vectorField>& motionPerPatch)")
        << "size of motion from bodyMotionFunction is not correct: " << motionPerPatch.size() << " != " << mesh_.boundaryMesh().size()
        << abort(FatalError);
    }

    //Loop over all patches in motionPerPatch from this body
    forAll(motionPerPatch,ipatch){
        //Ignore empty vectorfields
        if(motionPerPatch[ipatch].size() > 0){
            if(allMotionPerPatch[ipatch].size() == 0){
                allMotionPerPatch[ipatch] = vectorField(mesh_.Cf().boundaryField()[ipatch].size(),vector::zero);
            }
            allMotionPerPatch[ipatch] += motionPerPatch[ipatch];
        }
    }
}

void bodyCollector::setMovingPatchIDs()
{
    //Boolean list for tracking which patches are moving
    List<bool> patchIDsMoving(mesh_.boundaryMesh().size(),false);

    forAll(bodyMotions_,iBody){
        //Mark which patches are moving based on bodyMotionFunctions
        labelList movingPatchIDs = bmFunctions_[iBody]->motionPatchIDs();
        forAll(movingPatchIDs,ipatch){
            patchIDsMoving[movingPatchIDs[ipatch]] = true;
        }
    }

    //Get IDs of moving patches based on all bodyMotionFunctions
    wordList movingPatchNames(mesh_.boundaryMesh().size());
    label counter = 0;
    forAll(patchIDsMoving,ipatch){
        if(patchIDsMoving[ipatch]==true){
            movingPatchIDs_[counter] = ipatch;
            movingPatchNames[counter] = mesh_.boundaryMesh().names()[ipatch];
            counter++;
        }
    }
    movingPatchIDs_.setSize(counter);
    movingPatchNames.setSize(counter);

    Info << "bodyCollector found following moving patches: " << movingPatchNames << " with following patchIDs: " << movingPatchIDs_ << nl << endl;
}

// * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * //

bodyCollector::bodyCollector(const dictionary& dict,const fvMesh& mesh):
dict_(dict),
mesh_(mesh),
bmFunctions_(0),
bodyMotions_(dict.lookup("bodyMotions")),
movingPatchIDs_(mesh.boundaryMesh().size())
{
    //Boolean list to construct movingPatchIDs_ after all bodyMotionFunctions are checked for moving patches
    Info << nl << "Reading all bodyMotionFunctions and gathering information" << endl;

    bmFunctions_.setSize(bodyMotions_.size());
    forAll(bodyMotions_,iBody){
        //Name and type of current bodyMotion
        word bmName = bodyMotions_[iBody];
        word bmType = dict.subDict(bmName+"Dict").lookup("type");

        //Create autoPtr for this specific bodyMotionFunction
        autoPtr<bodyMotionFunction> bodyMotionFunctionPtr
        (
            bodyMotionFunction::New
            (
                bmType,
                bmName,
                dict,
                mesh
            )
        );

        bmFunctions_[iBody] = bodyMotionFunctionPtr;

    }

    setMovingPatchIDs();
}


// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> bodyCollector::getMotion()
{
    Info << nl << "Gathering motion from bodies" << endl;
    //Prepare return vector
    tmp<vectorField> tMotion(new vectorField(0));
    vectorField& allMotion = tMotion();

    //Create list to store all bodymotions in
    Field<vectorField> allMotionPerPatch(mesh_.boundaryMesh().size(),vectorField(0));

    //Loop over all bodyMotionFunctions to get their motion
    forAll(bmFunctions_,ibody){
        Field<vectorField> motionPerPatch = bmFunctions_[ibody]->getMotion();//Make this function in the body motion functions
        setAllMotion(allMotionPerPatch,motionPerPatch);
    }

    //Finally construct single vectorField with motion sorted in patchID order
    forAll(allMotionPerPatch,ipatch){
        //This should work, since all non moving patches should have zero size vectorFields
        allMotion.append(allMotionPerPatch[ipatch]);
    }

    Info << endl;

    return tMotion;
}

tmp<Field<vectorField> > bodyCollector::getPatchMotion()
{
    Info << nl << "Gathering motion from bodies" << endl;

    //Prepare return vector
    tmp<Field<vectorField> > tPatchMotion(new Field<vectorField>(mesh_.boundaryMesh().size(),vectorField(0)));
    Field<vectorField>& allMotionPerPatch = tPatchMotion();

    //Loop over all bodyMotionFunctions to get their motion
    forAll(bmFunctions_,ibody){
        Field<vectorField> motionPerPatch = bmFunctions_[ibody]->getMotion();//Make this function in the body motion functions
        setAllMotion(allMotionPerPatch,motionPerPatch);
    }

    Info << endl;
    return tPatchMotion;
}

tmp<Field<vectorField> > bodyCollector::getPatchMotion(const Field<vectorField>& patchForces)
{
    Info << nl << "Gathering motion from bodies based on fluid forces" << endl;
    //Prepare return vector
    tmp<Field<vectorField> > tPatchMotion(new Field<vectorField>(mesh_.boundaryMesh().size(),vectorField(0)));
    Field<vectorField>& allMotionPerPatch = tPatchMotion();

    //Loop over all bodyMotionFunctions to get their motion
    forAll(bmFunctions_,ibody){
        Field<vectorField> motionPerPatch = bmFunctions_[ibody]->getMotion(patchForces);//Make this function in the body motion functions
        setAllMotion(allMotionPerPatch,motionPerPatch);
    }
    Info << "Motion from bodies collected" << nl << endl;

    return tPatchMotion;
}

tmp<Field<vectorField> > bodyCollector::getSteadyStatePatchMotion(const Field<vectorField>& patchForces)
{
    Info << "bodyCollector::getSteadyStatePatchMotion(const Field<vectorField>& patchForces)" << endl;
    //Prepare return vector
    tmp<Field<vectorField> > tPatchMotion(new Field<vectorField>(mesh_.boundaryMesh().size(),vectorField(0)));
    Field<vectorField>& allMotionPerPatch = tPatchMotion();

    //Loop over all bodyMotionFunctions to get their motion
    forAll(bmFunctions_,ibody){
        Field<vectorField> motionPerPatch = bmFunctions_[ibody]->getSteadyStateMotion(patchForces);//Make this function in the body motion functions
        setAllMotion(allMotionPerPatch,motionPerPatch);
    }

    return tPatchMotion;
}

const labelList& bodyCollector::getMovingPatchIDs() const{
    return movingPatchIDs_;
}

void bodyCollector::update(){
    //Loop over all bodyMotionFunctions to call update functino
    forAll(bmFunctions_,ibody){
        bmFunctions_[ibody]->update();
    }
}

void bodyCollector::write(){
    forAll(bmFunctions_,ibody){
        bmFunctions_[ibody]->write();
    }
}

} // End namespace Foam

// ************************************************************************* //
