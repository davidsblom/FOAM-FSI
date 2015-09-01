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

#include "movingWallPressureGradient3FvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

Foam::scalarField Foam::movingWallPressureGradient3FvPatchScalarField::coeffs_( const volVectorField & vf ) const
{
    scalarField coeffs( 4, 0.0 );

    // if first time step
    if ( vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex() || vf.oldTime().oldTime().timeIndex() < 0 )
    {
        coeffs[0] = 1.0;
        coeffs[1] = -1.0;
    }
    else
    if ( vf.oldTime().oldTime().timeIndex() == vf.oldTime().oldTime().oldTime().timeIndex() || vf.oldTime().oldTime().oldTime().timeIndex() < 0 )
    {
        coeffs[0] = 3.0 / 2.0;
        coeffs[1] = -2.0;
        coeffs[2] = 1.0 / 2.0;
    }
    else
    {
        coeffs = coeffs_();
    }

    return coeffs;
}

Foam::scalarField Foam::movingWallPressureGradient3FvPatchScalarField::coeffs_() const
{
    scalarField coeffs( 4, 0.0 );
    coeffs[0] = 11.0 / 6.0;
    coeffs[1] = -3.0;
    coeffs[2] = 3.0 / 2.0;
    coeffs[3] = -1.0 / 3.0;

    return coeffs;
}

Foam::vectorField Foam::movingWallPressureGradient3FvPatchScalarField::ddtU(const volVectorField& U) const
{
    scalar deltaT = this->patch().boundaryMesh().mesh().time().deltaT().value();

    const scalarField coeffs = coeffs_( U );

    const scalar& patchID = this->patch().index();

    vectorField ddtU = (coeffs[0] * U.boundaryField()[patchID] + coeffs[1] * U.oldTime().boundaryField()[patchID] + coeffs[2] * U.oldTime().oldTime().boundaryField()[patchID] + coeffs[3] * U.oldTime().oldTime().oldTime().boundaryField()[patchID])/deltaT;

    return ddtU;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallPressureGradient3FvPatchScalarField::movingWallPressureGradient3FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho")
{}


Foam::movingWallPressureGradient3FvPatchScalarField::movingWallPressureGradient3FvPatchScalarField
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


Foam::movingWallPressureGradient3FvPatchScalarField::movingWallPressureGradient3FvPatchScalarField
(
    const movingWallPressureGradient3FvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::movingWallPressureGradient3FvPatchScalarField::movingWallPressureGradient3FvPatchScalarField
(
    const movingWallPressureGradient3FvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::movingWallPressureGradient3FvPatchScalarField::movingWallPressureGradient3FvPatchScalarField
(
    const movingWallPressureGradient3FvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallPressureGradient3FvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volVectorField& U =
        db().lookupObject<volVectorField>(UName_);

    //const volVectorField& gradPbc =
    //    db().lookupObject<volVectorField>("gradPbc");

    gradient() = -(ddtU(U) & patch().Sf())/patch().magSf();
    //const scalar& patchID = this->patch().index();

    //gradient() = (gradPbc.boundaryField()[patchID] & patch().Sf())/patch().magSf();

    //Info << "gradient = " << max(mag(gradient())) << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::movingWallPressureGradient3FvPatchScalarField::write(Ostream& os) const
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
        movingWallPressureGradient3FvPatchScalarField
    );
}

// ************************************************************************* //
