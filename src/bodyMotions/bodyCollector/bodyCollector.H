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
    bodyCollector.C
    newbodyCollector.C

\*---------------------------------------------------------------------------*/

#ifndef bodyCollector_H
#define bodyCollector_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"
#include "vectorField.H"
#include "foamTime.H"
#include "bodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bodyCollector Declaration
\*---------------------------------------------------------------------------*/

class bodyCollector
{
private:

    const dictionary& dict_;
    const fvMesh& mesh_;

    List<autoPtr<bodyMotionFunction> > bmFunctions_;

    wordList bodyMotions_;

    //List of patch IDs which are moving based on bodyMotinoFunctions
    labelList movingPatchIDs_;

    // Private Member Functions

        void setAllMotion(Field<vectorField>& allMotionPerPatch,const Field<vectorField>& motionPerPatch) const;

        //Set movingPatchIDs_ based on all bodyMotionFunctions
        void setMovingPatchIDs();

        //- Disallow copy construct
        bodyCollector(const bodyCollector&);

        //- Disallow default bitwise assignment
        void operator=(const bodyCollector&);

protected:


public:

    //- Runtime type information
    ClassName("bodyCollector");

    // Constructors
        bodyCollector(const dictionary& dict,const fvMesh& mesh);

    // Destructor
        ~bodyCollector()
        {}


    // Member Functions

        tmp<vectorField> getMotion();

        tmp<Field<vectorField> > getPatchMotion();

        tmp<Field<vectorField> > getPatchMotion(const Field<vectorField>& patchForces);

        tmp<Field<vectorField> > getSteadyStatePatchMotion(const Field<vectorField>& patchForces);

        const labelList& getMovingPatchIDs() const;

        void update();

        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
