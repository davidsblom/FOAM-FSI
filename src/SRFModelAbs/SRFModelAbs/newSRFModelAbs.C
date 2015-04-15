/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "SRFModelAbs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SRF
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<SRFModelAbs> SRFModelAbs::New
(
    const volVectorField& U
)
{
    word SRFModelAbsTypeName;

    // Enclose the creation of the SRFPropertiesDict to ensure it is
    // deleted before the SRFModelAbs is created - otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary SRFPropertiesDict
        (
            IOobject
            (
                "SRFProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        SRFPropertiesDict.lookup("SRFModelAbs") >> SRFModelAbsTypeName;
    }

    Info<< "Selecting SRFModelAbs " << SRFModelAbsTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(SRFModelAbsTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "SRFModelAbs::New(const fvMesh&)"
        )   << "Unknown SRFModelAbs type " << SRFModelAbsTypeName
            << nl << nl
            << "Valid SRFModelAbs types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<SRFModelAbs>(cstrIter()(U));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SRF
} // End namespace Foam

// ************************************************************************* //
