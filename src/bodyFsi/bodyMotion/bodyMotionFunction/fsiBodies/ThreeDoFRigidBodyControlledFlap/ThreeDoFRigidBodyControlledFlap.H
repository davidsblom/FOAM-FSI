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

Class
    ThreeDoFRigidBodyControlledFlap

Description
	Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    ThreeDoFRigidBodyControlledFlap.C

\*---------------------------------------------------------------------------*/

#ifndef ThreeDoFRigidBodyControlledFlap_H
#define ThreeDoFRigidBodyControlledFlap_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"
#include "Time.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ThreeDoFRigidBodyControlledFlap Declaration
\*---------------------------------------------------------------------------*/

class ThreeDoFRigidBodyControlledFlap
:
    public bodyMotionFunction
{
	// Private data
	const dictionary& dict_;
	const fvMesh& mesh_;
    const word name_;

    //Patches of this body
    wordList patchNames_;
    labelList patchIDs_;
    labelList flapPatches_;
    label flapPatchTop_;
    label flapPatchBottom_;

    //Solve properties
    word integrationScheme_;
    word flapIntegrationScheme_;

	//Body properties
	vector k_;
    vector c_;
    scalar ktheta_;
    scalar ctheta_;
    scalar m_;
    scalar Icg_;
    scalar Lcg2rc_;
    vector rc_;
    scalar thetaGeom_;//In radians

    //Flap properties
    vector topRotationOrigin_;
    vector bottomRotationOrigin_;
    scalar controllerStartTime_;
    scalar controllerVelocityGain_;
    scalar controllerAccGain_;
    
	//Body movement parameters: bodyState = [x,y,theta,dxdt,dydt,dthetadt]
    List<scalar> bodyState_;
    List<scalar> bodyStatePrevIter_;
    List<scalar> bodyStateOld_;
    List<scalar> bodyStateOldOld_;
    List<scalar> bodyState0_;
    List<scalar> bodyStateVel_;
    List<scalar> bodyStateVelOld_;
    scalar theta0_;
    scalar beta_;
    scalar betaOld_;
    scalar betaOldOld_;
    bool flapUpdated_;

    //Body forces
    vector bodyForces_;

    //Positions
    Field<vectorField> currentPosition_; 
    Field<vectorField> prevPosition_;
    Field<vectorField> initialPoints_;

    //Writing properties
    OFstream ofBody_;

    // Private Member Functions

        //- Writing function
        void writeState(const scalar& time);

        //- Calcaulte moment around current rotation center (rc_ + (x,y,rc_.z))
        vector calculateForcesAndMoments(const Field<vectorField>& patchForces) const;

        vector calculateWingPosition(const Field<vectorField>& patchForces);

        tmp<Field<vectorField> > calculateFlapMotion();

        //- Read initial state from dict if present (state = [x,y,theta,dxdt,dydt,dthetadt])
        void readInitialStateIfPresent();

        //- Set patchIDs_ based on the patchNames
        void setPatchIDs();

        //- Set rotation points based on shared point with non-flap patches
        void setFlapRotationPoints();

        //-
        point findRotationPoints(const labelList wingPatchIDs,const label flapPatchIndex) const;

    	//- Disallow default bitwise copy construct
        ThreeDoFRigidBodyControlledFlap(const ThreeDoFRigidBodyControlledFlap&);

    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        virtual tmp<Field<vectorField> > getChildMotion(const Field<vectorField>& patchForces);

        void setInitialPoints();

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;
	    
	public:

    //- Runtime type information
    TypeName("ThreeDoFRigidBodyControlledFlap");

	//Public data

    // Constructors

		//- Construct from dict + mesh + name
        ThreeDoFRigidBodyControlledFlap(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new ThreeDoFRigidBodyControlledFlap(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~ThreeDoFRigidBodyControlledFlap();

    // Member Functions

        //- Calculate steady state position based on forces
        //virtual vector calculateSteadyStatePosition(const vectorField& forces);

        virtual const Field<vectorField>& getPosition() const;

        virtual void update();

        //- Write out body properties if requested
        virtual void write();
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
