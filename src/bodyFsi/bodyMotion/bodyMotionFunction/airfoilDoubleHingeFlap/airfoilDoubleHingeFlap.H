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
    airfoilDoubleHingeFlap

Description
    Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    airfoilDoubleHingeFlap.C

\*---------------------------------------------------------------------------*/

#ifndef airfoilDoubleHingeFlap_H
#define airfoilDoubleHingeFlap_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dictionary.H"
#include "vector.H"
#include "vectorField.H"
#include "tmp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class airfoilDoubleHingeFlap Declaration
\*---------------------------------------------------------------------------*/

class airfoilDoubleHingeFlap
:
    public bodyMotionFunction
{
    private:

        //Useful references and name
        const dictionary& dict_;
        const fvMesh& mesh_;
        const word name_;

        //Patches of this body
        wordList flapPatches_;
        labelList flapPatchIDs_;
        label patchIndTop_;
        label patchIndBottom_;
        wordList wingPatches_;
        labelList wingPatchIDs_;
        wordList patchNames_;
        labelList patchIDs_;

        //Kinematic variables airfoil
    	scalar rotationAmplitude_;
    	scalar rotationFrequency_;
        vector rotationOrigin_;
        scalar rotationPhaseShift_;
    	scalar translationAmplitude_;
    	scalar translationFrequency_;
    	vector translationDirection_;
    	scalar translationPhaseShift_;

        //Kinematic variables flap
        scalar flapRotationAmplitude_;
        scalar flapRotationFrequency_;
        vector topRotationOrigin_;
        vector bottomRotationOrigin_;
        labelList topRotationPointIDs_;
        labelList bottomRotationPointIDs_;
        scalar flapPhaseShift_;
        scalar transitionFraction_;

        //General kinematic settings
        bool smoothStart_;

        //Positions
        Field<vectorField> currentPosition_;
        Field<vectorField> prevPosition_;
        Field<vectorField> initialPoints_;

        //Writing file
        OFstream ofBody_;

        //- Disallow default bitwise copy construct
        airfoilDoubleHingeFlap(const airfoilDoubleHingeFlap&);

        //- Writing function
        void writeState(const scalar& time);

        //- Calculate flap angle and returns it
        scalar getFlapAngle(const scalar time) const;

        //- Calculate wing translation (x,y) and rotation angle and return in vector (x,y,angle)
        vector getWingState(const scalar& time) const;

        //Private functions
        tmp<Field<vectorField> > calculatePosition(const scalar time);

        //- Set rotation points based on shared point with non-flap patches
        void setFlapRotationPoints();

        //- Find rotation points based on the top and bottom flap patches
        labelList findRotationPoints(const labelList wingPatchIDs,const label flapPatchIndex) const;

        //Set patch IDs of flap
        void setFlapPatchIDs();

        //Set patch IDs of wing
        void setWingPatchIDs();

        //Set patch IDs of flap and wing in patchNames_ and patchIDs_
        void setPatchNamesAndIDs();


    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        void setInitialPoints();

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;

    public:

    //- Runtime type information
    TypeName("airfoilDoubleHingeFlap");

    //Public data

    // Constructors

        //- Construct from dict + mesh
        airfoilDoubleHingeFlap(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new airfoilDoubleHingeFlap(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~airfoilDoubleHingeFlap();

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
