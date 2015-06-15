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
    dummyFsiBody

Description
    Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    dummyFsiBody.C

\*---------------------------------------------------------------------------*/

#ifndef dummyFsiBody_H
#define dummyFsiBody_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"
#include "Time.H"
#include "simpleMatrix.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dummyFsiBody Declaration
\*---------------------------------------------------------------------------*/

class dummyFsiBody
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

    Field<vectorField> initialPoints_;

    //Writing properties
    OFstream ofBody_;

    // Private Member Functions

        //- Writing function
        void writeState(const scalar& time);

        
        void setInitialPoints();

        //- Set patchIDs_ based on the patchNames
        void setPatchIDs();

        //- Disallow default bitwise copy construct
        dummyFsiBody(const dummyFsiBody&);

    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        virtual tmp<Field<vectorField> > getChildMotion(const Field<vectorField>& patchForces);

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;

    public:

    //- Runtime type information
    TypeName("dummyFsiBody");

    //Public data

    // Constructors

        //- Construct from dict + mesh + name
        dummyFsiBody(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new dummyFsiBody(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~dummyFsiBody();

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
