/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kkLOmegaof24.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kkLOmegaof24, 0);
addToRunTimeSelectionTable(RASModel, kkLOmegaof24, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> kkLOmegaof24::fv(const volScalarField& Ret) const
{
    return(1.0 - exp(-sqrt(Ret)/Av_));
}


tmp<volScalarField> kkLOmegaof24::fINT() const
{
    return
    (
        min
        (
            kt_/(Cint_*(kl_ + kt_ + k0_)),
            dimensionedScalar("1.0", dimless, 1.0)
        )
    );
}


tmp<volScalarField> kkLOmegaof24::fSS(const volScalarField& Omega) const
{
    return(exp(-sqr(Css_*nu()*Omega/(kt_ + k0_))));
}


tmp<volScalarField> kkLOmegaof24::Cmu(const volScalarField& S) const
{
    return(1.0/(A0_ + As_*(S/(omega_ + omega0_))));
}


tmp<volScalarField> kkLOmegaof24::BetaTS(const volScalarField& ReOmega) const
{
    return(scalar(1) - exp(-sqr(max(ReOmega - CtsCrit_, scalar(0)))/Ats_));
}


tmp<volScalarField> kkLOmegaof24::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -CtauL_*ktL
          /
            (
                sqr
                (
                    lambdaEff*Omega
                  + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        dimLength*inv(dimTime),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}


tmp<volScalarField> kkLOmegaof24::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return(fv*CmuStd_*sqrt(ktS)*lambdaEff);
}


tmp<volScalarField> kkLOmegaof24::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        scalar(1)
      - exp
        (
           -0.41
           *pow4
            (
                lambdaEff
              / (
                    lambdaT
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        lambdaT.dimensions(),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}


tmp<volScalarField> kkLOmegaof24::phiBP(const volScalarField& Omega) const
{
    return
    (
        min
        (
            max
            (
                kt_/nu()
             / (
                    Omega
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        Omega.dimensions(),
                        ROOTVSMALL
                    )
                )
              - CbpCrit_,
                scalar(0)
            ),
            scalar(50.0)
        )
    );
}


tmp<volScalarField> kkLOmegaof24::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return
    (
        max
        (
            ReOmega
          - CnatCrit_
          / (
                fNatCrit + dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    );
}


tmp<volScalarField> kkLOmegaof24::D(const volScalarField& k) const
{
    return nu()*magSqr(fvc::grad(sqrt(k)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kkLOmegaof24::kkLOmegaof24
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
            coeffDict_,
            2.12
        )
    ),
    Av_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Av",
            coeffDict_,
            6.75
        )
    ),
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            coeffDict_,
            1250
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            coeffDict_,
            1000
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
            coeffDict_,
            0.02
        )
    ),
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            coeffDict_,
            3.4e-6
        )
    ),
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            coeffDict_,
            0.12
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            coeffDict_,
            0.035
        )
    ),
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
            coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.92
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            0.3
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            coeffDict_,
            1.5
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            coeffDict_,
            2.495
        )
    ),
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
            coeffDict_,
            0.09
        )
    ),
    Prtheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prtheta",
            coeffDict_,
            0.85
        )
    ),
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
            coeffDict_,
            1
        )
    ),
    Sigmaw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaw",
            coeffDict_,
            1.17
        )
    ),
    kt_
    (
        IOobject
        (
            "kt",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    kl_
    (
        IOobject
        (
            "kl",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_
        ),
        kt_*omega_ + D(kl_) + D(kt_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    y_(mesh_)
{
    bound(kt_, k0_);
    bound(kl_, k0_);
    bound(omega_, omega0_);
    bound(epsilon_, epsilon0_);

    // Evaluating nut_ is complex so start from the field read from file
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kkLOmegaof24::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k() - nut_*twoSymm(fvc::grad(U_)),
            kt_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kkLOmegaof24::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kkLOmegaof24::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      //- fvc::div(nuEff()*dev(T(fvc::grad(U))))
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

bool kkLOmegaof24::read()
{
    if (RASModel::read())
    {
        A0_.readIfPresent(coeffDict());
        As_.readIfPresent(coeffDict());
        Av_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Anat_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Ats_.readIfPresent(coeffDict());
        CbpCrit_.readIfPresent(coeffDict());
        Cnc_.readIfPresent(coeffDict());
        CnatCrit_.readIfPresent(coeffDict());
        Cint_.readIfPresent(coeffDict());
        CtsCrit_.readIfPresent(coeffDict());
        CrNat_.readIfPresent(coeffDict());
        C11_.readIfPresent(coeffDict());
        C12_.readIfPresent(coeffDict());
        CR_.readIfPresent(coeffDict());
        CalphaTheta_.readIfPresent(coeffDict());
        Css_.readIfPresent(coeffDict());
        CtauL_.readIfPresent(coeffDict());
        Cw1_.readIfPresent(coeffDict());
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        CwR_.readIfPresent(coeffDict());
        Clambda_.readIfPresent(coeffDict());
        CmuStd_.readIfPresent(coeffDict());
        Prtheta_.readIfPresent(coeffDict());
        Sigmak_.readIfPresent(coeffDict());
        Sigmaw_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kkLOmegaof24::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
        y_.boundaryField() = max(y_.boundaryField(), VSMALL);
    }

    const volScalarField lambdaT(sqrt(kt_)/(omega_ + omega0_));

    const volScalarField lambdaEff(min(Clambda_*y_, lambdaT));

    const volScalarField fw
    (
        pow
        (
            lambdaEff
           /(lambdaT + dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2.0/3.0
        )
    );

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));

    const volScalarField S2(2.0*magSqr(dev(symm(gradU))));

    const volScalarField ktS(fSS(Omega)*fw*kt_);

    const volScalarField nuts
    (
        fv(sqr(fw)*kt_/nu()/(omega_ + omega0_))
       *fINT()
       *Cmu(sqrt(S2))*sqrt(ktS)*lambdaEff
    );
    const volScalarField Pkt(nuts*S2);

    const volScalarField ktL(kt_ - ktS);
    const volScalarField ReOmega(sqr(y_)*Omega/nu());
    const volScalarField nutl
    (
        min
        (
            C11_*fTaul(lambdaEff, ktL, Omega)*Omega*sqr(lambdaEff)
           *sqrt(ktL)*lambdaEff/nu()
          + C12_*BetaTS(ReOmega)*ReOmega*sqr(y_)*Omega
        ,
            0.5*(kl_ + ktL)/(sqrt(S2) + omega0_)
        )
    );

    const volScalarField Pkl(nutl*S2);

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff, fv(sqr(fw)*kt_/nu()/(omega_ + omega0_)), ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("SMALL", dimless, ROOTVSMALL);

    const volScalarField Rbp
    (
        CR_*(1.0 - exp(-phiBP(Omega)()/Abp_))*omega_
       /(fw + fwMin)
    );

    const volScalarField fNatCrit(1.0 - exp(-Cnc_*sqrt(kl_)*y_/nu()));

    // Natural source term divided by kl_
    const volScalarField Rnat
    (
        CrNat_*(1.0 - exp(-phiNAT(ReOmega, fNatCrit)/Anat_))*Omega
    );


    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(alphaTEff), omega_)
     ==
        Cw1_*Pkt*omega_/(kt_ + k0_)
      - fvm::SuSp
        (
            (1.0 - CwR_/(fw + fwMin))*kl_*(Rbp + Rnat)/(kt_ + k0_)
          , omega_
        )
      - fvm::Sp(Cw2_*sqr(fw)*omega_, omega_)
      + (
            Cw3_*fOmega(lambdaEff, lambdaT)*alphaTEff*sqr(fw)*sqrt(kt_)
        )().dimensionedInternalField()/pow3(y_.dimensionedInternalField())
    );

    omegaEqn().relax();
    //omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omega0_);


    const volScalarField Dl(D(kl_));

    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(kl_)
      + fvm::div(phi_, kl_)
      - fvm::laplacian(nu(), kl_)
     ==
        Pkl
      - fvm::Sp(Rbp + Rnat + Dl/(kl_ + k0_), kl_)
    );

    klEqn().relax();
    //klEqn().boundaryManipulate(kl_.boundaryField());

    solve(klEqn);
    bound(kl_, k0_);


    const volScalarField Dt(D(kt_));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(kt_)
      + fvm::div(phi_, kt_)
      - fvm::laplacian(DkEff(alphaTEff), kt_)
     ==
        Pkt
      + (Rbp + Rnat)*kl_
      - fvm::Sp(omega_ + Dt/(kt_+ k0_), kt_)
    );

    ktEqn().relax();
    //ktEqn().boundaryManipulate(kt_.boundaryField());

    solve(ktEqn);
    bound(kt_, k0_);


    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = kt_*omega_ + Dl + Dt;
    bound(epsilon_, epsilon0_);


    // Re-calculate turbulent viscosity
    nut_ = nuts + nutl;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
