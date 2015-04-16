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

#include "SRFModelAbs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace SRF
    {
        defineTypeNameAndDebug(SRFModelAbs, 0);
        defineRunTimeSelectionTable(SRFModelAbs, dictionary);
    }
}


void Foam::SRF::SRFModelAbs::setPatchProperties()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelHashSet excludedPatches(excludedPatchLabels_);

    excludedPatches_ = labelList(patches.size(), 0);
    includedPatches_ = labelList(patches.size(), 0);
    label nIncludedPatches=0;
    label nExcludedPatches=0;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.isWall() && !excludedPatches.found(patchI))
        {
            //Add to included patches
            includedPatches_[nIncludedPatches]=patchI;
            nIncludedPatches++;
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            //Add to excluded patches
            excludedPatches_[nExcludedPatches]=patchI;
            nExcludedPatches++;
        }
    }
    includedPatches_.setSize(nIncludedPatches);
    excludedPatches_.setSize(nExcludedPatches);

    //Info << nl <<  "Patch selection for SRF complete:" << endl;
    //Info << "includedPatches = " << includedPatches_ << endl;
    //Info << "excludedPatches = " << excludedPatches_ << nl << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SRF::SRFModelAbs::SRFModelAbs
(
    const word& type,
    const volVectorField& U
)
:
    IOdictionary
    (
        IOobject
        (
            "SRFProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    U_(U),
    mesh_(U_.mesh()),
    axis_(lookup("axis")),
    origin_(lookup("origin")),
    SRFModelAbsCoeffs_(subDict(type + "Coeffs")),
    omega_(dimensionedVector("omega", dimless/dimTime, vector::zero)),
    excludedPatchNames_
    (
        lookupOrDefault("nonRotatingPatches", wordList(0))
    )
{
    // Normalise the axis
    if (mag(axis_) < SMALL)
    {
        FatalErrorIn
        (
            "SRF::SRFModelAbs::SRFModelAbs\n"
            "(\n"
            "    const word& type,\n"
            "    const volVectorField& Urel\n"
            ")"
        )   << "Zero length axis: " << axis_ << ".  This is not allowed."
            << abort(FatalError);
    }

    axis_ /= mag(axis_);

    //Get ID from excluded patches
    excludedPatchLabels_.setSize(excludedPatchNames_.size());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(excludedPatchNames_, i)
    {
        excludedPatchLabels_[i] = patches.findPatchID(excludedPatchNames_[i]);

        if (excludedPatchLabels_[i] == -1)
        {
            FatalErrorIn
            (
            "Foam::SRF::SRFModelAbs::SRFModelAbs(const word&, const volVectorField&)"
            )   << "cannot find SRF patch " << excludedPatchNames_[i] << " from nonRotatingPatches"
            << exit(FatalError);
        }
    }

    setPatchProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SRF::SRFModelAbs::~SRFModelAbs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SRF::SRFModelAbs::read()
{
    if (regIOobject::read())
    {
        // Re-read axis
        lookup("axis") >> axis_;
        lookup("origin") >> origin_;

        if (mag(axis_) < SMALL)
        {
            FatalErrorIn("SRF::SRFModelAbs::read()")
                << "Zero length axis: " << axis_ << ".  This is not allowed."
                << abort(FatalError);
        }

        axis_ /= mag(axis_);

        // Re-read sub-model coeffs
        SRFModelAbsCoeffs_ = subDict(type() + "Coeffs");

        //Get ID from excluded patches
        excludedPatchLabels_.setSize(excludedPatchNames_.size());
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        forAll(excludedPatchNames_, i)
        {
            excludedPatchLabels_[i] = patches.findPatchID(excludedPatchNames_[i]);

            if (excludedPatchLabels_[i] == -1)
            {
                FatalErrorIn
                (
                "Foam::SRF::SRFModelAbs::SRFModelAbs(const word&, const volVectorField&)"
                )   << "cannot find SRF patch " << excludedPatchNames_[i] << " from nonRotatingPatches"
                << exit(FatalError);
            }
        }

        setPatchProperties();

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::vector& Foam::SRF::SRFModelAbs::axis() const
{
    return axis_;
}


const Foam::dimensionedVector& Foam::SRF::SRFModelAbs::omega() const
{
    return omega_;
}

void Foam::SRF::SRFModelAbs::correctBoundaryVelocity(volVectorField& U) const
{
    Info << "SRFModelAbs::correctBoundaryVelocity(volVectorField& U)" << endl;
    //Included patches
    const vector& omega = omega_.value();

    forAll(includedPatches_, patchi)
    {
        const label& patchID = includedPatches_[patchi];
        const vectorField& patchC = mesh_.Cf().boundaryField()[patchID];
        //const vectorField patchR = (patchC - axis_*(axis_ & patchC));
        const vectorField patchR = (patchC - origin_);

        U.boundaryField()[patchID] == (omega ^ patchR);
    }
}

void Foam::SRF::SRFModelAbs::addCoriolis(fvVectorMatrix& UEqn) const
{
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    const vector& omega = omega_.value();

    forAll(V, i)
    {
        Usource[i] -= V[i]*(omega ^ U[i]);
    }
}



Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh> >
Foam::SRF::SRFModelAbs::Fcoriolis() const
{
    return tmp<DimensionedField<vector, volMesh> >
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                "Fcoriolis",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (omega_ ^ U_)
        )
    );
}

void Foam::SRF::SRFModelAbs::relativeFlux(surfaceScalarField& phi) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();
    const vector& omega = omega_.value();

    //phi.internalField() -= (omega ^ (Cf.internalField() - axis_*(axis_ & Cf.internalField()))) & Sf.internalField();
    phi.internalField() -= (omega ^ (Cf.internalField() - origin_)) & Sf.internalField();

    // Included patches
    forAll(includedPatches_, patchi)
    {
        const label& patchID = includedPatches_[patchi];
        phi.boundaryField()[patchID] = 0.0;
    }

    // Excluded patches
    forAll(excludedPatches_, patchi)
    {
        const label& patchID = excludedPatches_[patchi];
        //phi.boundaryField()[patchID] -=
        //    (omega ^ (Cf.boundaryField()[patchID] - axis_*(axis_ & Cf.boundaryField()[patchID]))) & Sf.boundaryField()[patchID];
        phi.boundaryField()[patchID] -= (omega ^ (Cf.boundaryField()[patchID] - origin_)) & Sf.boundaryField()[patchID];
    }
}

void Foam::SRF::SRFModelAbs::absoluteFlux(surfaceScalarField& phi) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();
    const vector& omega = omega_.value();

    //phi.internalField() += (omega ^ (Cf.internalField() - axis_*(axis_ & Cf.internalField()))) & Sf.internalField();
    phi.internalField() += (omega ^ (Cf.internalField() - origin_)) & Sf.internalField();

    // Included patches
    forAll(includedPatches_, patchi)
    {
        const label& patchID = includedPatches_[patchi];
        //phi.boundaryField()[patchID] +=
        //    (omega ^ (Cf.boundaryField()[patchID] - axis_*(axis_ & Cf.boundaryField()[patchID]))) & Sf.boundaryField()[patchID];
        phi.boundaryField()[patchID] += (omega ^ (Cf.boundaryField()[patchID] - origin_)) & Sf.boundaryField()[patchID];
    }

    // Excluded patches
    forAll(excludedPatches_, patchi)
    {
        const label& patchID = excludedPatches_[patchi];
        //phi.boundaryField()[patchID] +=
        //    (omega ^ (Cf.boundaryField()[patchID] - axis_*(axis_ & Cf.boundaryField()[patchID]))) & Sf.boundaryField()[patchID];
        phi.boundaryField()[patchID] += (omega ^ (Cf.boundaryField()[patchID] - origin_)) & Sf.boundaryField()[patchID];
    }
}
// ************************************************************************* //
