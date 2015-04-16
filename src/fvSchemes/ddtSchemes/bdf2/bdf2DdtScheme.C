
/*---------------------------------------------------------------------------*\
*  =========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\    /   O peration     |
\\  /    A nd           | Copyright held by original author
\\\\\\\\\\/     M anipulation  |
\\ -------------------------------------------------------------------------------
\\ License
\\  This file is part of OpenFOAM.
\\
\\  OpenFOAM is free software; you can redistribute it and/or modify it
\\  under the terms of the GNU General Public License as published by the
\\  Free Software Foundation; either version 2 of the License, or (at your
\\  option) any later version.
\\
\\  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
\\  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
\\  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
\\  for more details.
\\
\\  You should have received a copy of the GNU General Public License
\\  along with OpenFOAM; if not, write to the Free Software Foundation,
\\  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\\
\*---------------------------------------------------------------------------*/

#include "bdf2DdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "basicSymmetryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace fv
    {
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        template<class Type>
        scalar bdf2DdtScheme<Type>::deltaT_() const
        {
            return mesh().time().deltaT().value();
        }

        template<class Type>
        scalar bdf2DdtScheme<Type>::deltaT0_() const
        {
            return mesh().time().deltaT0().value();
        }

        template<class Type>
        template<class GeoField>
        scalar bdf2DdtScheme<Type>::deltaT0_( const GeoField & vf ) const
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

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf2DdtScheme<Type>::fvcDdt
            ( const dimensioned<Type> & dt )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddt(" + dt.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_();

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            if ( mesh().moving() )
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
                            dt.dimensions() / dimTime,
                            pTraits<Type>::zero
                        )
                    )
                );

                tdtdt().internalField() = rDeltaT.value() * dt.value() *
                    (
                    coefft - ( coefft0 * mesh().V0() - coefft00 * mesh().V00() ) / mesh().V()
                    );

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
                            dt.dimensions() / dimTime,
                            pTraits<Type>::zero
                        ),
                        calculatedFvPatchField<Type>::typeName
                    )
                       );
            }
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf2DdtScheme<Type>::fvcDdt
            ( const GeometricField<Type, fvPatchField, volMesh> & vf )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddt(" + vf.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            if ( mesh().moving() )
            {
                return tmp<GeometricField<Type, fvPatchField, volMesh> >
                       (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        ddtIOobject,
                        mesh(),
                        rDeltaT.dimensions() * vf.dimensions(),
                        rDeltaT.value() *
                        (
                            coefft * vf.internalField() -
                            (
                                coefft0 * vf.oldTime().internalField() * mesh().V0()
                                - coefft00 * vf.oldTime().oldTime().internalField()
                                * mesh().V00()
                            ) / mesh().V()
                        ),
                        rDeltaT.value() *
                        (
                            coefft * vf.boundaryField() -
                            (
                                coefft0 * vf.oldTime().boundaryField()
                                - coefft00 * vf.oldTime().oldTime().boundaryField()
                            )
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
                        rDeltaT *
                        (
                            coefft * vf
                            - coefft0 * vf.oldTime()
                            + coefft00 * vf.oldTime().oldTime()
                        )
                    )
                       );
            }
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf2DdtScheme<Type>::fvcDdt
        (
            const dimensionedScalar & rho,
            const GeometricField<Type, fvPatchField, volMesh> & vf
        )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddt(" + rho.name() + ',' + vf.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            if ( mesh().moving() )
            {
                return tmp<GeometricField<Type, fvPatchField, volMesh> >
                       (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        ddtIOobject,
                        mesh(),
                        rDeltaT.dimensions() * rho.dimensions() * vf.dimensions(),
                        rDeltaT.value() * rho.value() *
                        (
                            coefft * vf.internalField() -
                            (
                                coefft0 * vf.oldTime().internalField() * mesh().V0()
                                - coefft00 * vf.oldTime().oldTime().internalField()
                                * mesh().V00()
                            ) / mesh().V()
                        ),
                        rDeltaT.value() * rho.value() *
                        (
                            coefft * vf.boundaryField() -
                            (
                                coefft0 * vf.oldTime().boundaryField()
                                - coefft00 * vf.oldTime().oldTime().boundaryField()
                            )
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
                        rDeltaT * rho *
                        (
                            coefft * vf
                            - coefft0 * vf.oldTime()
                            + coefft00 * vf.oldTime().oldTime()
                        )
                    )
                       );
            }
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf2DdtScheme<Type>::fvcDdt
        (
            const volScalarField & rho,
            const GeometricField<Type, fvPatchField, volMesh> & vf
        )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddt(" + rho.name() + ',' + vf.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            if ( mesh().moving() )
            {
                return tmp<GeometricField<Type, fvPatchField, volMesh> >
                       (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        ddtIOobject,
                        mesh(),
                        rDeltaT.dimensions() * rho.dimensions() * vf.dimensions(),
                        rDeltaT.value() *
                        (
                            coefft * rho.internalField() * vf.internalField() -
                            (
                                coefft0 * rho.oldTime().internalField()
                                * vf.oldTime().internalField() * mesh().V0()
                                - coefft00 * rho.oldTime().oldTime().internalField()
                                * vf.oldTime().oldTime().internalField() * mesh().V00()
                            ) / mesh().V()
                        ),
                        rDeltaT.value() *
                        (
                            coefft * rho.boundaryField() * vf.boundaryField() -
                            (
                                coefft0 * rho.oldTime().boundaryField()
                                * vf.oldTime().boundaryField()
                                - coefft00 * rho.oldTime().oldTime().boundaryField()
                                * vf.oldTime().oldTime().boundaryField()
                            )
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
                        rDeltaT *
                        (
                            coefft * rho * vf
                            - coefft0 * rho.oldTime() * vf.oldTime()
                            + coefft00 * rho.oldTime().oldTime() * vf.oldTime().oldTime()
                        )
                    )
                       );
            }
        }

        template<class Type>
        tmp<fvMatrix<Type> >
        bdf2DdtScheme<Type>::fvmDdt
            ( GeometricField<Type, fvPatchField, volMesh> & vf )
        {
            tmp<fvMatrix<Type> > tfvm
            (
                new fvMatrix<Type>
                (
                    vf,
                    vf.dimensions() * dimVol / dimTime
                )
            );

            fvMatrix<Type> & fvm = tfvm();

            scalar rDeltaT = 1.0 / deltaT_();

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            fvm.diag() = (coefft * rDeltaT) * mesh().V();

            if ( mesh().moving() )
            {
                fvm.source() = rDeltaT *
                    (
                    coefft0 * vf.oldTime().internalField() * mesh().V0()
                    - coefft00 * vf.oldTime().oldTime().internalField()
                    * mesh().V00()
                    );
            }
            else
            {
                fvm.source() = rDeltaT * mesh().V() *
                    (
                    coefft0 * vf.oldTime().internalField()
                    - coefft00 * vf.oldTime().oldTime().internalField()
                    );
            }

            return tfvm;
        }

        template<class Type>
        tmp<fvMatrix<Type> >
        bdf2DdtScheme<Type>::fvmDdt
        (
            const dimensionedScalar & rho,
            GeometricField<Type, fvPatchField, volMesh> & vf
        )
        {
            tmp<fvMatrix<Type> > tfvm
            (
                new fvMatrix<Type>
                (
                    vf,
                    rho.dimensions() * vf.dimensions() * dimVol / dimTime
                )
            );
            fvMatrix<Type> & fvm = tfvm();

            scalar rDeltaT = 1.0 / deltaT_();

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            fvm.diag() = ( coefft * rDeltaT * rho.value() ) * mesh().V();

            if ( mesh().moving() )
            {
                fvm.source() = rDeltaT * rho.value() *
                    (
                    coefft0 * vf.oldTime().internalField() * mesh().V0()
                    - coefft00 * vf.oldTime().oldTime().internalField()
                    * mesh().V00()
                    );
            }
            else
            {
                fvm.source() = rDeltaT * mesh().V() * rho.value() *
                    (
                    coefft0 * vf.oldTime().internalField()
                    - coefft00 * vf.oldTime().oldTime().internalField()
                    );
            }

            return tfvm;
        }

        template<class Type>
        tmp<fvMatrix<Type> >
        bdf2DdtScheme<Type>::fvmDdt
        (
            const volScalarField & rho,
            GeometricField<Type, fvPatchField, volMesh> & vf
        )
        {
            tmp<fvMatrix<Type> > tfvm
            (
                new fvMatrix<Type>
                (
                    vf,
                    rho.dimensions() * vf.dimensions() * dimVol / dimTime
                )
            );
            fvMatrix<Type> & fvm = tfvm();

            scalar rDeltaT = 1.0 / deltaT_();

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            fvm.diag() = (coefft * rDeltaT) * rho.internalField() * mesh().V();

            if ( mesh().moving() )
            {
                fvm.source() = rDeltaT *
                    (
                    coefft0 * rho.oldTime().internalField()
                    * vf.oldTime().internalField() * mesh().V0()
                    - coefft00 * rho.oldTime().oldTime().internalField()
                    * vf.oldTime().oldTime().internalField() * mesh().V00()
                    );
            }
            else
            {
                fvm.source() = rDeltaT * mesh().V() *
                    (
                    coefft0 * rho.oldTime().internalField()
                    * vf.oldTime().internalField()
                    - coefft00 * rho.oldTime().oldTime().internalField()
                    * vf.oldTime().oldTime().internalField()
                    );
            }

            return tfvm;
        }

        template<class Type>
        tmp<typename bdf2DdtScheme<Type>::fluxFieldType>
        bdf2DdtScheme<Type>::fvcDdtPhiCorr
        (
            const volScalarField & rA,
            const GeometricField<Type, fvPatchField, volMesh> & U,
            const fluxFieldType & phi
        )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( U );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

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
                dimensioned<scalar>( "1", dimless, 1.0 )
            );

            forAll( U.boundaryField(), patchI )
            {
                ddtPhiCoeff.boundaryField()[patchI] = 0.0;
            }

            // === make sure fvc::interpolate is performed on 1.0/rA seperately === //
            return tmp<fluxFieldType>
                   (
                new fluxFieldType
                (
                    ddtIOobject,
                    rDeltaT * ddtPhiCoeff // fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
                    * (
                        (
                            (
                                coefft0 * phi.oldTime()
                                - coefft00 * phi.oldTime().oldTime()
                            )
                            - (
                                fvc::interpolate
                                (
                                    (
                                        coefft0 * U.oldTime()
                                        - coefft00 * U.oldTime().oldTime()
                                    )
                                ) & mesh().Sf()
                                )
                        ) / fvc::interpolate( 1.0 / rA )
                        )
                )
                   );
        }

        // Specialized member function for ddtPhiCorr when type = vector
        template<>
        tmp<surfaceScalarField> bdf2DdtScheme<vector>::fvcDdtPhiCorr
        (
            const volScalarField & rA,
            const GeometricField<vector, fvPatchField, volMesh> & U,
            const surfaceScalarField & phi
        )
        {
            // Info << "Consistent backwardDdtPhiCorr" << endl;

            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( U );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0); // = 3/2
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) ); // = 1/2
            scalar coefft0 = coefft + coefft00; // = 2, which should be -2 for good sign in scheme itself

            // Set coefficient for excluding certain boundaries from ddtPhiCorr
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
                dimensioned<scalar>( "1", dimless, 1.0 )
            );

            // Set coefficients at boundary which have fixedValues to 0
            forAll( U.boundaryField(), patchI )
            {
                if
                (
                    U.boundaryField()[patchI].fixesValue()

                    /*|| isA<symmetryFvPatchVectorField>(U.boundaryField()[patchI])
                    || isA<basicSymmetryFvPatchVectorField>(U.boundaryField()[patchI])
                    || isA<slipFvPatchVectorField>(U.boundaryField()[patchI])*/
                )
                {
                    ddtPhiCoeff.boundaryField()[patchI] = 0.0;
                }
            }

            // If mesh is moving
            if ( mesh().moving() )
            {
                // Set volume ratio of oldV over V
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

                V0oV.internalField() = mesh().V0() / mesh().V();
                V0oV.correctBoundaryConditions();

                // Set volume ratio of oldoldV over V
                volScalarField V00oV
                (
                    IOobject
                    (
                        "V00oV",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimless,
                    zeroGradientFvPatchScalarField::typeName
                );

                V00oV.internalField() = mesh().V00() / mesh().V();
                V00oV.correctBoundaryConditions();

                // Dimension check
                if
                (
                    U.dimensions() == dimVelocity
                    && phi.dimensions() == dimVelocity * dimArea
                )
                {
                    // Get Uf with old times recorded
                    const surfaceVectorField & Uf = mesh().objectRegistry::lookupObject<surfaceVectorField>( "Uf" );

                    // Construct parts of H which need to be subtracted
                    surfaceVectorField UV0oV = -fvc::interpolate( V0oV * U.oldTime() ); // term coming from H/A
                    surfaceVectorField UV00oV = -fvc::interpolate( V00oV * U.oldTime().oldTime() ); // term coming from H/A

                    // Create return variable
                    tmp<surfaceScalarField> ddtPhiCorr
                    (
                        new surfaceScalarField
                        (
                            ddtIOobject,
                            rDeltaT * ddtPhiCoeff
                            * (
                                coefft0 * ( ( UV0oV + fvc::interpolate( V0oV ) * ( Uf.oldTime() ) ) & mesh().Sf() )
                                - coefft00 * ( ( UV00oV + fvc::interpolate( V00oV ) * ( Uf.oldTime().oldTime() ) ) & mesh().Sf() )
                                ) / fvc::interpolate( 1.0 / rA )
                        )
                    );

                    return ddtPhiCorr;
                }
                else
                {
                    FatalErrorIn
                    (
                        "backwardDdtScheme<vector>::fvcDdtPhiCorr"
                    ) << "dimensions of phi are not correct"
                      << abort( FatalError );

                    return fluxFieldType::null();
                }
            }
            else
            {
                // If not moving simplify calculations
                if
                (
                    U.dimensions() == dimVelocity
                    && phi.dimensions() == dimVelocity * dimArea
                )
                {
                    surfaceScalarField phi0 = phi.oldTime() - ( fvc::interpolate( U.oldTime() ) & mesh().Sf() );
                    surfaceScalarField phi00 = phi.oldTime().oldTime() - ( fvc::interpolate( U.oldTime().oldTime() ) & mesh().Sf() );

                    // Create return variable
                    return tmp<fluxFieldType>
                           (
                        new fluxFieldType
                        (
                            ddtIOobject,
                            rDeltaT * ddtPhiCoeff
                            * (
                                coefft0 * phi0
                                - coefft00 * phi00
                                )
                            / fvc::interpolate( 1.0 / rA )
                        )
                           );
                }
                else
                {
                    FatalErrorIn
                    (
                        "backwardDdtScheme<vector>::fvcDdtPhiCorr"
                    ) << "dimensions of phi are not correct"
                      << abort( FatalError );

                    return fluxFieldType::null();
                }
            }
        }

        template<class Type>
        tmp<typename bdf2DdtScheme<Type>::fluxFieldType>
        bdf2DdtScheme<Type>::fvcDdtPhiCorr
        (
            const volScalarField & rA,
            const volScalarField & rho,
            const GeometricField<Type, fvPatchField, volMesh> & U,
            const fluxFieldType & phiAbs
        )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

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

            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( U );

            scalar coefft = 1 + deltaT / (deltaT + deltaT0);
            scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
            scalar coefft0 = coefft + coefft00;

            if
            (
                U.dimensions() == dimVelocity
                && phiAbs.dimensions() == dimVelocity * dimArea
            )
            {
                return tmp<fluxFieldType>
                       (
                    new fluxFieldType
                    (
                        ddtIOobject,
                        rDeltaT * this->fvcDdtPhiCoeff( U.oldTime(), phiAbs.oldTime() )
                        * (
                            coefft0 * fvc::interpolate( rA * rho.oldTime() )
                            * phiAbs.oldTime()
                            - coefft00 * fvc::interpolate( rA * rho.oldTime().oldTime() )
                            * phiAbs.oldTime().oldTime()
                            - (
                                fvc::interpolate
                                (
                                    rA *
                                    (
                                        coefft0 * rho.oldTime() * U.oldTime()
                                        - coefft00 * rho.oldTime().oldTime()
                                        * U.oldTime().oldTime()
                                    )
                                ) & mesh().Sf()
                                )
                            )
                    )
                       );
            }
            else
            if
            (
                U.dimensions() == dimVelocity
                && phiAbs.dimensions() == rho.dimensions() * dimVelocity * dimArea
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
                            phiAbs.oldTime() / fvc::interpolate( rho.oldTime() )
                        )
                        * (
                            fvc::interpolate( rA * rho.oldTime() )
                            * (
                                coefft0 * phiAbs.oldTime()
                                / fvc::interpolate( rho.oldTime() )
                                - coefft00 * phiAbs.oldTime().oldTime()
                                / fvc::interpolate( rho.oldTime().oldTime() )
                                )
                            - (
                                fvc::interpolate
                                (
                                    rA * rho.oldTime() *
                                    (
                                        coefft0 * U.oldTime()
                                        - coefft00 * U.oldTime().oldTime()
                                    )
                                ) & mesh().Sf()
                                )
                            )
                    )
                       );
            }
            else
            if
            (
                U.dimensions() == rho.dimensions() * dimVelocity
                && phiAbs.dimensions() == rho.dimensions() * dimVelocity * dimArea
            )
            {
                return tmp<fluxFieldType>
                       (
                    new fluxFieldType
                    (
                        ddtIOobject,
                        rDeltaT
                        * this->fvcDdtPhiCoeff( rho.oldTime(), U.oldTime(), phiAbs.oldTime() )
                        * (
                            fvc::interpolate( rA )
                            * (
                                coefft0 * phiAbs.oldTime()
                                - coefft00 * phiAbs.oldTime().oldTime()
                                )
                            - (
                                fvc::interpolate
                                (
                                    rA *
                                    (
                                        coefft0 * U.oldTime()
                                        - coefft00 * U.oldTime().oldTime()
                                    )
                                ) & mesh().Sf()
                                )
                            )
                    )
                       );
            }
            else
            {
                FatalErrorIn
                (
                    "bdf2DdtScheme<Type>::fvcDdtPhiCorr"
                ) << "dimensions of phiAbs are not correct"
                  << abort( FatalError );

                return fluxFieldType::null();
            }
        }

        template<class Type>
        tmp<surfaceScalarField> bdf2DdtScheme<Type>::meshPhi
            ( const GeometricField<Type, fvPatchField, volMesh> & vf )
        {
            scalar deltaT = deltaT_();
            scalar deltaT0 = deltaT0_( vf );

            // Coefficient for t-3/2 (between times 0 and 00)
            scalar coefft0_00 = deltaT / (deltaT + deltaT0);

            // Coefficient for t-1/2 (between times n and 0)
            scalar coefftn_0 = 1 + coefft0_00;

            return coefftn_0 * mesh().phi() - coefft0_00 * mesh().phi().oldTime();
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace fv

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
