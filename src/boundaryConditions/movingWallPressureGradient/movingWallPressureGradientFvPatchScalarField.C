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

#include "movingWallPressureGradientFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

Foam::scalar Foam::movingWallPressureGradientFvPatchScalarField::deltaT_() const
{
    return this->patch().boundaryMesh().mesh().time().deltaT().value();
}

Foam::scalar Foam::movingWallPressureGradientFvPatchScalarField::deltaT0_() const
{
    return this->patch().boundaryMesh().mesh().time().deltaT0().value();
}

Foam::scalar Foam::movingWallPressureGradientFvPatchScalarField::deltaT0_( const volVectorField& vf ) const
{
    // Bug fix, Zeljko Tukovic: solver with outer iterations over a time-step
    // HJ, 12/Feb/2010
    // Bug fix, Thijs Gillebaart: in first time step -> solve backwardEuler insteadt of BDF2. Especially for meshPhi.
    // Because mesh.phi.oldtime is zero when not yet defined. Constant velocity causes different wall velocity at first time step otherwise
    // Solved by using vf.oldTime().oldTime().timeIndex() < 0
    // TG, 22/Feb/2013
    if ( vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex() || vf.oldTime().oldTime().timeIndex() < 0 )
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}

Foam::vectorField Foam::movingWallPressureGradientFvPatchScalarField::ddtU(const volVectorField& U) const
{
    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_( U );

    scalar coefft = 1 + deltaT / (deltaT + deltaT0);
    scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
    scalar coefft0 = coefft + coefft00;

    const scalar& patchID = this->patch().index();

    vectorField ddtU = coefft * U.boundaryField()[patchID] - coefft0 * U.oldTime().boundaryField()[patchID] + coefft00 * U.oldTime().oldTime().boundaryField()[patchID];
    return ddtU;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallPressureGradientFvPatchScalarField::movingWallPressureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho")
{}


Foam::movingWallPressureGradientFvPatchScalarField::movingWallPressureGradientFvPatchScalarField
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


Foam::movingWallPressureGradientFvPatchScalarField::movingWallPressureGradientFvPatchScalarField
(
    const movingWallPressureGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::movingWallPressureGradientFvPatchScalarField::movingWallPressureGradientFvPatchScalarField
(
    const movingWallPressureGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


Foam::movingWallPressureGradientFvPatchScalarField::movingWallPressureGradientFvPatchScalarField
(
    const movingWallPressureGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallPressureGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volVectorField& U =
        db().lookupObject<volVectorField>(UName_);

    const fvPatchField<vector>& Up =
        lookupPatchField<volVectorField, vector>(UName_);

    /*if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        phip /= rhop;
    }*/


    gradient() = (ddtU(U) & patch().Sf())/patch().magSf();
    //gradient() = (phip - (patch().Sf() & Up))/patch().magSf()/rAp;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::movingWallPressureGradientFvPatchScalarField::write(Ostream& os) const
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
        movingWallPressureGradientFvPatchScalarField
    );
}

// ************************************************************************* //
