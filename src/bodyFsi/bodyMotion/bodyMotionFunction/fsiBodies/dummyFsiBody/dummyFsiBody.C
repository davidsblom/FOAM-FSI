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
    dummyFsiBody

\*---------------------------------------------------------------------------*/
#include "dummyFsiBody.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(dummyFsiBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, dummyFsiBody, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
void dummyFsiBody::writeState(const scalar& time)
{
    //Body
    ofBody_ << time << "\t" << 0 << endl;
}

void dummyFsiBody::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
    }
}

void dummyFsiBody::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void dummyFsiBody::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

// * * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

tmp<Field<vectorField> > dummyFsiBody::getChildMotion()
{
    FatalErrorIn("void dummyFsiBody::getChildMotion()")
            << "this bodyMotionFunction cannot be used in a non-fsi solver, since getChildMotion(const Field<vectorField>& patchForces) needs to be called"
            << abort(FatalError);

    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));

    return tMotion;
}

tmp<Field<vectorField> > dummyFsiBody::getChildMotion(const Field<vectorField>& patchForces)
{
    //Prepare return variable
    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    forAll(patchIDs_,ipatch){
        motion[ipatch] = vectorField(initialPoints_[ipatch].size(),vector::zero);
    }

    return tMotion;
}

const labelList dummyFsiBody::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList dummyFsiBody::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dummyFsiBody::dummyFsiBody(const dictionary& dict,const fvMesh& mesh,const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    setPatchIDs();
    setInitialPoints();

    if(writeToFile())
    {
        ofBody_ << "time\tbody" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dummyFsiBody::~dummyFsiBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& dummyFsiBody::getPosition() const{
    return initialPoints_;
}

void dummyFsiBody::update()
{
    bodyMotionFunction::update();//This updates localTimeIndex
}

void dummyFsiBody::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
