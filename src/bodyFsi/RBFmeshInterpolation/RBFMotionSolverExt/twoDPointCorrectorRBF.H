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
    Foam::twoDPointCorrectorRBF

Description
    Class applies a two-dimensional correction to mesh motion point field.

    The correction guarantees that the mesh does not get twisted during motion
    and thus introduce a third dimension into a 2-D problem.

    The operation is performed by looping through all edges approximately
    normal to the plane and enforcing their orthoginality onto the plane by
    adjusting points on their ends.

SourceFiles
    twoDPointCorrectorRBF.C

\*---------------------------------------------------------------------------*/

#ifndef twoDPointCorrectorRBF_H
#define twoDPointCorrectorRBF_H

#include "pointField.H"
#include "labelList.H"
#include "vector.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class twoDPointCorrectorRBF Declaration
\*---------------------------------------------------------------------------*/

class twoDPointCorrectorRBF
:
    public twoDPointCorrector
{
    // Private data

        //- Reference to moving mesh
        const polyMesh& mesh_;

        //- Holder for marker value: -1 if on "back" side and 0 otherwise
        labelList pointMarker_;

        //- Point IDs marked with 0
        labelList useablePointIDs_;

        //- Point IDs marked with -1
        labelList shadowPointIDs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        twoDPointCorrectorRBF(const twoDPointCorrectorRBF&);

        //- Disallow default bitwise assignment
        void operator=(const twoDPointCorrectorRBF&);

        void setMarker();

public:

    // Constructors

        //- Construct from components
        twoDPointCorrectorRBF(const polyMesh& mesh);


    // Destructor

        ~twoDPointCorrectorRBF();


    // Member Functions

        const labelList& marker() const;

        void setShadowSide(vectorField& newpoints) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
