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
    twodRigidBody

Description
	Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    twodRigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef twodRigidBody_H
#define twodRigidBody_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dictionary.H"
#include "vector.H"
#include "vectorField.H"
#include "tmp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twodRigidBody Declaration
\*---------------------------------------------------------------------------*/

class twodRigidBody
:
    public bodyMotionFunction
{
    private:

    //Useful references and name
	const dictionary& dict_;
	const fvMesh& mesh_;
    const word name_;

    //Patches of this body
    wordList patchNames_;
    labelList patchIDs_;

    //Kinematic variables
	scalar rotationAmplitude_;
	scalar rotationFrequency_;
    vector rotationOrigin_;
    scalar rotationPhaseShift_;
	scalar translationAmplitude_;
	scalar translationFrequency_;
	vector translationDirection_;
	scalar translationPhaseShift_;
    bool smoothStart_;

    //Positions
    Field<vectorField> currentPosition_;
    Field<vectorField> prevPosition_;
    Field<vectorField> initialPoints_;

    //Writing properties
    OFstream ofBody_;

        // Private Member Functions

        //- Writing function
        void writeState(const scalar& time);

        //- Get state back for time given (x,y,angle)
        vector getState(const scalar& time);

    	//- Disallow default bitwise copy construct
        twodRigidBody(const twodRigidBody&);

        //Private functions
        tmp<Field<vectorField> > calculatePosition(const scalar time);

        void setPatchIDs();

    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        void setInitialPoints();

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;

	public:

    //- Runtime type information
    TypeName("twodRigidBody");

	//Public data

    // Constructors

        //- Construct from dict + mesh
		twodRigidBody(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new twodRigidBody(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~twodRigidBody();

    // Member Functions

		virtual const Field<vectorField>& getPosition() const;

        virtual void update();

        //- Write out body properties if requested
        virtual void write();
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
