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
    RBFFunctionTG

Description
    Radial basis function virtual base class

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    RBFFunctionTG.C
    newRBFFunctionTG.C

\*---------------------------------------------------------------------------*/

#ifndef RBFFunctionTG_H
#define RBFFunctionTG_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RBFFunctionTG Declaration
\*---------------------------------------------------------------------------*/

class RBFFunctionTG
{
    // Private Member Functions

        //- Disallow copy construct
        RBFFunctionTG(const RBFFunctionTG&);

        //- Disallow default bitwise assignment
        void operator=(const RBFFunctionTG&);


public:

    //- Runtime type information
    TypeName("RBFFunctionTG");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RBFFunctionTG,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Return a pointer to the selected RBF function
        static autoPtr<RBFFunctionTG> New
        (
            const word& type,
            const dictionary& dict
        );


    // Constructors

        //- Construct null
        RBFFunctionTG()
        {}

        //- Create and return a clone
        virtual autoPtr<RBFFunctionTG> clone() const = 0;


    // Destructor

        virtual ~RBFFunctionTG()
        {}


    // Member Functions

        //- Return RBF weights
        virtual tmp<scalarField> weights
        (
            const vectorField& points,
            const vector& controlPoint
        ) const = 0;
        
        //- Return RBF weights for control matrix
        virtual tmp<scalarField> weights
        (
        	const vectorField& points
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
