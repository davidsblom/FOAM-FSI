/*---------------------------------------------------------------------------* \
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

\*---------------------------------------------------------------------------*/

#include "EulerRCDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerRCDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt().internalField() =
            rDeltaT.value()*dt.value()*(1.0 - mesh().V0()/mesh().V());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerRCDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerRCDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerRCDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
EulerRCDdtScheme<Type>::fvmDdt
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerRCDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerRCDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<typename EulerRCDdtScheme<Type>::fluxFieldType>
EulerRCDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phiAbs.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

	// === Set boundaries correct of U === //
    surfaceScalarField ddtPhiCoeff
    (
        IOobject
        (
            "ddtPhiCoeff",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensioned<scalar>("1", dimless, 1.0)
    );

    forAll (U.boundaryField(), patchI)
    {
        ddtPhiCoeff.boundaryField()[patchI] = 0.0;
    }
	
    //tmp<fluxFieldType> phiCorr =
    //    phiAbs.oldTime() - (fvc::interpolate(U.oldTime()) & mesh().Sf());
    tmp<fluxFieldType> phiCorr =
        phiAbs.oldTime() - (fvc::interpolate(U.oldTime()) & mesh().Sf());//TEMP: Everything should be old including Sf (Paper tukovic about Rhie-Chow interpolation)

	// === make sure fvc::interpolate is performed on 1.0/rA seperately === //
    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            rDeltaT * ddtPhiCoeff//fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime(), phiCorr())
           *1.0/fvc::interpolate(1.0/rA)*phiCorr
        )
    );
}

template<>
tmp<surfaceScalarField> EulerRCDdtScheme<vector>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const surfaceScalarField& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    //Set coefficient for excluding certain boundaries from ddtPhiCorr
    surfaceScalarField ddtPhiCoeff
    (
        IOobject
        (
            "ddtPhiCoeff",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensioned<scalar>("1", dimless, 1.0)
    );

    //Set coefficients at boundary which have fixedValues to 0
    forAll(U.boundaryField(), patchI)
    {
        if
        (
            U.boundaryField()[patchI].fixesValue()
         //|| isA<symmetryFvPatchVectorField>(U.boundaryField()[patchI])
         //|| isA<slipFvPatchVectorField>(U.boundaryField()[patchI])
        )
        {
            ddtPhiCoeff.boundaryField()[patchI] = 0.0;
        }
    }

    //If mesh is moving
    if (mesh().moving())
    {
        //Set volume ratio of oldV over V
        volScalarField V0oV
        (
            IOobject
            (
                "V0oV",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V0oV.internalField() = mesh().V0()/mesh().V();
        V0oV.correctBoundaryConditions();

        //Dimension check
        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            //Set old velocity interpolated to new mesh
            surfaceVectorField U0 = fvc::interpolate(U.oldTime());

            //Get Sf with old times recored
            const surfaceVectorField& Sf = 
                mesh().objectRegistry::lookupObject<surfaceVectorField>
                (
                    "Sf"
                );

            //Set ddtPhiCorr step for Uold
            U0 += Sf.oldTime()
               *(phi.oldTime() - (Sf.oldTime()&U0))
               /(
                    magSqr(Sf.oldTime())
                  + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
                );
            
            //Create return variable
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        (
                            (
                                fvc::interpolate(V0oV)*U0
                            ) & mesh().Sf()
                        ) 
                      - (
                            fvc::interpolate
                            (
                                V0oV*U.oldTime()
                            ) & mesh().Sf()
                        )
                    )/fvc::interpolate(1.0/rA)
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "EulerDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
    else
    {
        //If not moving simplify calculations
        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            //Set interpolated old velocity
            surfaceVectorField dU0 =  fvc::interpolate(U.oldTime());
            forAll(dU0.boundaryField(), patchI)
            {
                if (!U.boundaryField()[patchI].coupled())
                {
                    dU0.boundaryField()[patchI] = 
                        U.oldTime().boundaryField()[patchI]
                       .patchInternalField();
                }
            }
            
            //Set ddtPhiCorr (difference between flux and Uold & Sf) for t(n)
            dU0 = mesh().Sf()
               *(phi.oldTime() - (mesh().Sf()&dU0))
               /(
                    magSqr(mesh().Sf())
                  + dimensionedScalar("SMALL", sqr(dimArea), SMALL)
                );

            //Create return variable
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    rDeltaT*ddtPhiCoeff
                   *(
                        mesh().Sf()&dU0
                    )
                   /fvc::interpolate(1.0/rA) 
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "EulerDdtScheme<vector>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<typename EulerRCDdtScheme<Type>::fluxFieldType>
EulerRCDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ','
      + rho.name() + ','
      + U.name() + ','
      + phiAbs.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               * this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime())
               *(
                   fvc::interpolate(rA*rho.oldTime())*phiAbs.oldTime()
                 - (fvc::interpolate(rA*rho.oldTime()*U.oldTime())
                  & mesh().Sf())
                )
            )
        );
    }
    else if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               * this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                )
               *(
                   fvc::interpolate(rA*rho.oldTime())
                  *phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                 - (
                       fvc::interpolate
                       (
                           rA*rho.oldTime()*U.oldTime()
                       ) & mesh().Sf()
                   )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               * this->fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phiAbs.oldTime())
               *(
                   fvc::interpolate(rA)*phiAbs.oldTime()
                 - (fvc::interpolate(rA*U.oldTime()) & mesh().Sf())
                )
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "EulerRCDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of phiAbs are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> EulerRCDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return mesh().phi();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
