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
    fileTranslationalBody

Description
	Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    fileTranslationalBody.C

\*---------------------------------------------------------------------------*/

#ifndef fileTranslationalBody_H
#define fileTranslationalBody_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dictionary.H"
#include "vector.H"
#include "vectorField.H"
#include "tmp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fileTranslationalBody Declaration
\*---------------------------------------------------------------------------*/

class fileTranslationalBody
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
    string fileName_;
	scalar scaleFactor_;
	scalar frequency_;
	vector translationDirection_;
    scalarField amplitudes_;
    scalar deltaT_;
    label timeIndex_;
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
        scalar getState(const scalar& time);

    	//- Disallow default bitwise copy construct
        fileTranslationalBody(const fileTranslationalBody&);

        //Private functions
        tmp<Field<vectorField> > calculatePosition(const scalar time);

        void setPatchIDs();

        //Read in file
        scalarField readFile(string fileName) const;

        void checkTimeStepSettings() const;

    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        void setInitialPoints();

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;

	public:

    //- Runtime type information
    TypeName("fileTranslationalBody");

	//Public data

    // Constructors

        //- Construct from dict + mesh
		fileTranslationalBody(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new fileTranslationalBody(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~fileTranslationalBody();

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
