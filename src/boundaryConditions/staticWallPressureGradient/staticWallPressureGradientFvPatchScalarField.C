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

\*---------------------------------------------------------------------------*/

#include "staticWallPressureGradientFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

Foam::vectorField Foam::staticWallPressureGradientFvPatchScalarField::getLapNuUpatch(const volVectorField& lapNuU) const
{
    const scalar& patchID = this->patch().index();

    vectorField lapNuUpatch = lapNuU.boundaryField()[patchID];

    return lapNuUpatch;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticWallPressureGradientFvPatchScalarField::staticWallPressureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho")
{}


Foam::staticWallPressureGradientFvPatchScalarField::staticWallPressureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::staticWallPressureGradientFvPatchScalarField::staticWallPressureGradientFvPatchScalarField
(
    const staticWallPressureGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::staticWallPressureGradientFvPatchScalarField::staticWallPressureGradientFvPatchScalarField
(
    const staticWallPressureGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::staticWallPressureGradientFvPatchScalarField::staticWallPressureGradientFvPatchScalarField
(
    const staticWallPressureGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::staticWallPressureGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volVectorField& lapNuU = db().lookupObject<volVectorField>("lapNuU");

    //const fvPatchField<vector>& Up = lookupPatchField<volVectorField, vector>(UName_);
    
    gradient() = -(getLapNuUpatch(lapNuU) & patch().Sf())/patch().magSf();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::staticWallPressureGradientFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        staticWallPressureGradientFvPatchScalarField
    );
}

// ************************************************************************* //
