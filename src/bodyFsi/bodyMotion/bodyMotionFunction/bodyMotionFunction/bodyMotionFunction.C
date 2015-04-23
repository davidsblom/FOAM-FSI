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
bodyMotionFunction

\*---------------------------------------------------------------------------*/

#include "bodyMotionFunction.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(bodyMotionFunction, 0);
    defineRunTimeSelectionTable(bodyMotionFunction, dictionary);

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

	void bodyMotionFunction::sortChildMotion(Field<vectorField>& motion,const labelList& childPatchIDs,const Field<vectorField>& childMotion) const{
		forAll(childMotion,iChildPatches){
			label patchID = childPatchIDs[iChildPatches];//patchID of this specific entry
			if(patchID >= 0 && patchID < mesh_.boundaryMesh().size()){
				motion[patchID] = childMotion[iChildPatches];//set motion in list at patchID location
			}else{
				FatalErrorIn("tmp<Field<vectorField> > bodyMotionFunction::sortChildMotion()")
			    << "patchID from child is not in range: " << patchID << " (0:" << mesh_.boundaryMesh().size() << ")"
			    << abort(FatalError);
			}
		}
	}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //

	const label& bodyMotionFunction::localTimeIndex() const{
		return timeIndex_;
	}

	tmp<Field<vectorField> > bodyMotionFunction::getChildMotion(const Field<vectorField>& forces){
		return getChildMotion();
	}

	tmp<Field<vectorField> > bodyMotionFunction::getChildSteadyStateMotion(const Field<vectorField>& forces){
		return getChildMotion(forces);
	}

// * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * //

	bodyMotionFunction::bodyMotionFunction(const dictionary& dict, const fvMesh& mesh,const word name):
	name_(name),
	mesh_(mesh),
	timeIndex_(-1),
	write_(dict.lookupOrDefault("write",false))
	{

	}

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

    tmp<Field<vectorField> > bodyMotionFunction::getMotion()
    {
    	//Prepare tmp return value
		tmp<Field<vectorField> > tMotion(new Field<vectorField>(mesh_.boundaryMesh().size(),vectorField(0)));
		Field<vectorField>& motion = tMotion();

		//Get the motion and the order in which they are sorted in the List<vectorField>
		Field<vectorField> childMotion = getChildMotion();
		labelList childPatchIDs = getChildMotionPatchIDs();

		//Set childMotion in appropiate place of motion
		sortChildMotion(motion,childPatchIDs,childMotion);

		return tMotion;
    }

    tmp<Field<vectorField> > bodyMotionFunction::getMotion(const Field<vectorField>& patchForces)
    {
    	//Prepare tmp return value
		tmp<Field<vectorField> > tMotion(new Field<vectorField>(mesh_.boundaryMesh().size(),vectorField(0)));
		Field<vectorField>& motion = tMotion();

		//Get the motion and the order in which they are sorted in the List<vectorField>
		Field<vectorField> childMotion = getChildMotion(patchForces);
		labelList childPatchIDs = getChildMotionPatchIDs();

		//Set childMotion in appropiate place of motion
		sortChildMotion(motion,childPatchIDs,childMotion);

		return tMotion;
    }

    tmp<Field<vectorField> > bodyMotionFunction::getSteadyStateMotion(const Field<vectorField>& patchForces)
    {
    	//Prepare tmp return value
		tmp<Field<vectorField> > tMotion(new Field<vectorField>(mesh_.boundaryMesh().size(),vectorField(0)));
		Field<vectorField>& motion = tMotion();

		//Get the motion and the order in which they are sorted in the List<vectorField>
		Field<vectorField> childMotion = getChildSteadyStateMotion(patchForces);
		labelList childPatchIDs = getChildMotionPatchIDs();

		//Set childMotion in appropiate place of motion
		sortChildMotion(motion,childPatchIDs,childMotion);

		return tMotion;
    }

	const labelList bodyMotionFunction::motionPatchIDs() const{
		return getChildMotionPatchIDs();
	}

	const wordList bodyMotionFunction::motionPatchNames() const{
		return getChildMotionPatchNames();
	}

    word bodyMotionFunction::name() const{
    	return name_;
    }

    void bodyMotionFunction::update(){
    	timeIndex_ = mesh_.time().timeIndex();
    }

	void bodyMotionFunction::write(){

	}

    bool bodyMotionFunction::writeToFile() const{
    	return write_;
    }

} // End namespace Foam

// ************************************************************************* //
