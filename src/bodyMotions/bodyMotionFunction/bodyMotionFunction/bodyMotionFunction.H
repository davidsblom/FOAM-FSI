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
    rigidMotion

Description
    Rigid motion virtual base class

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    bodyMotionFunction.C
    newbodyMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef bodyMotionFunction_H
#define bodyMotionFunction_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "primitiveFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class bodyMotionFunction
{
private:

    //Name of body
    const word name_;

    //Refernce of mesh
    const fvMesh& mesh_;

    //Record time index to only return motion when timeIndex < mesh.timeIndex
    label timeIndex_;

    bool write_;

    // Private Member Functions

        void sortChildMotion(Field<vectorField>& motion,const labelList& childPatchIDs,const Field<vectorField>& childMotion) const;

        //- Disallow copy construct
        bodyMotionFunction(const bodyMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const bodyMotionFunction&);

protected:
        virtual tmp<Field<vectorField> > getChildMotion() = 0;

        virtual tmp<Field<vectorField> > getChildMotion(const Field<vectorField>& patchForces);

        virtual tmp<Field<vectorField> > getChildSteadyStateMotion(const Field<vectorField>& patchForces);

        virtual const labelList getChildMotionPatchIDs() const = 0;

        virtual const wordList getChildMotionPatchNames() const = 0;

        const label& localTimeIndex() const;

public:

    //- Runtime type information
    TypeName("bodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            bodyMotionFunction,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word name
            ),
            (dict, mesh, name)
        );


    // Selectors

        //- Return a pointer to the selected RigidMotion function
        static autoPtr<bodyMotionFunction> New
        (
            const word& type,
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct null
        //bodyMotionFunction()
        //{}

        bodyMotionFunction(const dictionary& dict, const fvMesh& mesh,const word name);

        //- Create and return a clone
        virtual autoPtr<bodyMotionFunction> clone() const = 0;


    // Destructor

        virtual ~bodyMotionFunction()
        {}


    // Member Functions

        //- Set position of wing based on time
        virtual const Field<vectorField>& getPosition() const = 0;

        tmp<Field<vectorField> > getMotion();

        tmp<Field<vectorField> > getMotion(const Field<vectorField>& patchForces);

        tmp<Field<vectorField> > getSteadyStateMotion(const Field<vectorField>& patchForces);

        const labelList motionPatchIDs() const;

        const wordList motionPatchNames() const;

        word name() const;

        virtual void update();

        //- Write out body properties if requested - empty function here. Can be filled by childeren if wanted
        virtual void write();

        //- Return boolean on writing
        bool writeToFile() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
