
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

#include "bdf3DdtScheme.H"
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
        scalar bdf3DdtScheme<Type>::deltaT_() const
        {
            return mesh().time().deltaT().value();
        }

        template<class Type>
        scalar bdf3DdtScheme<Type>::deltaT0_() const
        {
            return mesh().time().deltaT0().value();
        }

        template<class Type>
        template<class GeoField>
        scalar bdf3DdtScheme<Type>::deltaT0_( const GeoField & vf ) const
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

        template<class Type>
        template<class GeoField>
        scalarField bdf3DdtScheme<Type>::coeffs_( const GeoField & vf ) const
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

        template<class Type>
        scalarField bdf3DdtScheme<Type>::coeffs_() const
        {
            scalarField coeffs( 4, 0.0 );
            coeffs[0] = 11.0 / 6.0;
            coeffs[1] = -3.0;
            coeffs[2] = 3.0 / 2.0;
            coeffs[3] = -1.0 / 3.0;

            return coeffs;
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf3DdtScheme<Type>::fvcDdt
            ( const dimensioned<Type> & dt )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddt(" + dt.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalarField coeffs = coeffs_();

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

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
                    coeffs[0] +
                    (
                        coeffs[1] * mesh().V0()
                        + coeffs[2] * mesh().V00()
                        + coeffs[3] * V000.dimensionedInternalField()
                    ) / mesh().V()
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
        bdf3DdtScheme<Type>::fvcDdt
            ( const GeometricField<Type, fvPatchField, volMesh> & vf )
        {
            dimensionedScalar rDeltaT = 1.0 / mesh().time().deltaT();

            IOobject ddtIOobject
            (
                "ddt(" + vf.name() + ')',
                mesh().time().timeName(),
                mesh()
            );

            scalarField coeffs = coeffs_( vf );

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

                return tmp<GeometricField<Type, fvPatchField, volMesh> >
                       (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        ddtIOobject,
                        mesh(),
                        rDeltaT.dimensions() * vf.dimensions(),
                        rDeltaT.value() *
                        (
                            coeffs[0] * vf.internalField() +
                            (
                                coeffs[1] * vf.oldTime().internalField() * mesh().V0()
                                + coeffs[2] * vf.oldTime().oldTime().internalField() * mesh().V00()
                                + coeffs[3] * vf.oldTime().oldTime().oldTime().internalField() * V000.dimensionedInternalField()
                            ) / mesh().V()
                        ),
                        rDeltaT.value() *
                        (
                            coeffs[0] * vf.boundaryField()
                            + coeffs[1] * vf.oldTime().boundaryField()
                            + coeffs[2] * vf.oldTime().oldTime().boundaryField()
                            + coeffs[3] * vf.oldTime().oldTime().oldTime().boundaryField()
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
                            coeffs[0] * vf
                            + coeffs[1] * vf.oldTime()
                            + coeffs[2] * vf.oldTime().oldTime()
                            + coeffs[3] * vf.oldTime().oldTime().oldTime()
                        )
                    )
                       );
            }
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf3DdtScheme<Type>::fvcDdt
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

            scalarField coeffs = coeffs_( vf );

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

                return tmp<GeometricField<Type, fvPatchField, volMesh> >
                       (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        ddtIOobject,
                        mesh(),
                        rDeltaT.dimensions() * rho.dimensions() * vf.dimensions(),
                        rDeltaT.value() * rho.value() *
                        (
                            coeffs[0] * vf.internalField() +
                            (
                                coeffs[1] * vf.oldTime().internalField() * mesh().V0()
                                + coeffs[2] * vf.oldTime().oldTime().internalField() * mesh().V00()
                                + coeffs[3] * vf.oldTime().oldTime().oldTime().internalField() * V000.dimensionedInternalField()
                            ) / mesh().V()
                        ),
                        rDeltaT.value() * rho.value() *
                        (
                            coeffs[0] * vf.boundaryField()
                            + coeffs[1] * vf.oldTime().boundaryField()
                            + coeffs[2] * vf.oldTime().oldTime().boundaryField()
                            + coeffs[3] * vf.oldTime().oldTime().oldTime().boundaryField()
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
                            coeffs[0] * vf
                            + coeffs[1] * vf.oldTime()
                            + coeffs[2] * vf.oldTime().oldTime()
                            + coeffs[3] * vf.oldTime().oldTime().oldTime()
                        )
                    )
                       );
            }
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        bdf3DdtScheme<Type>::fvcDdt
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

            scalarField coeffs = coeffs_( vf );

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

                return tmp<GeometricField<Type, fvPatchField, volMesh> >
                       (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        ddtIOobject,
                        mesh(),
                        rDeltaT.dimensions() * rho.dimensions() * vf.dimensions(),
                        rDeltaT.value() *
                        (
                            coeffs[0] * rho.internalField() * vf.internalField() +
                            (
                                coeffs[1] * rho.oldTime().internalField() * vf.oldTime().internalField() * mesh().V0()
                                + coeffs[2] * rho.oldTime().oldTime().internalField() * vf.oldTime().oldTime().internalField() * mesh().V00()
                                + coeffs[3] * rho.oldTime().oldTime().oldTime().internalField() * vf.oldTime().oldTime().oldTime().internalField() * V000.dimensionedInternalField()
                            ) / mesh().V()
                        ),
                        rDeltaT.value() *
                        (
                            coeffs[0] * rho.boundaryField() * vf.boundaryField()
                            + coeffs[1] * rho.oldTime().boundaryField() * vf.oldTime().boundaryField()
                            + coeffs[2] * rho.oldTime().oldTime().boundaryField() * vf.oldTime().oldTime().boundaryField()
                            + coeffs[3] * rho.oldTime().oldTime().oldTime().boundaryField() * vf.oldTime().oldTime().oldTime().boundaryField()
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
                            coeffs[0] * rho * vf
                            + coeffs[1] * rho.oldTime() * vf.oldTime()
                            + coeffs[2] * rho.oldTime().oldTime() * vf.oldTime().oldTime()
                            + coeffs[3] * rho.oldTime().oldTime().oldTime() * vf.oldTime().oldTime().oldTime()
                        )
                    )
                       );
            }
        }

        template<class Type>
        tmp<fvMatrix<Type> >
        bdf3DdtScheme<Type>::fvmDdt
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

            scalarField coeffs = coeffs_( vf );

            fvm.diag() = (coeffs[0] * rDeltaT) * mesh().V();

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

                fvm.source() = rDeltaT *
                    (
                    -coeffs[1] * vf.oldTime().internalField() * mesh().V0()
                    - coeffs[2] * vf.oldTime().oldTime().internalField() * mesh().V00()
                    - coeffs[3] * vf.oldTime().oldTime().oldTime().internalField() * V000.dimensionedInternalField()
                    );
            }
            else
            {
                fvm.source() = rDeltaT * mesh().V() *
                    (
                    -coeffs[1] * vf.oldTime().internalField()
                    - coeffs[2] * vf.oldTime().oldTime().internalField()
                    - coeffs[3] * vf.oldTime().oldTime().oldTime().internalField()
                    );
            }

            return tfvm;
        }

        template<class Type>
        tmp<fvMatrix<Type> >
        bdf3DdtScheme<Type>::fvmDdt
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

            scalarField coeffs = coeffs_( vf );

            fvm.diag() = ( coeffs[0] * rDeltaT * rho.value() ) * mesh().V();

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

                fvm.source() = rDeltaT * rho.value() *
                    (
                    -coeffs[1] * vf.oldTime().internalField() * mesh().V0()
                    - coeffs[2] * vf.oldTime().oldTime().internalField() * mesh().V00()
                    - coeffs[3] * vf.oldTime().oldTime().oldTime().internalField() * V000.dimensionedInternalField()
                    );
            }
            else
            {
                fvm.source() = rDeltaT * mesh().V() * rho.value() *
                    (
                    -coeffs[1] * vf.oldTime().internalField()
                    - coeffs[2] * vf.oldTime().oldTime().internalField()
                    - coeffs[3] * vf.oldTime().oldTime().oldTime().internalField()
                    );
            }

            return tfvm;
        }

        template<class Type>
        tmp<fvMatrix<Type> >
        bdf3DdtScheme<Type>::fvmDdt
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

            scalarField coeffs = coeffs_( vf );

            fvm.diag() = (coeffs[0] * rDeltaT) * rho.internalField() * mesh().V();

            if ( mesh().moving() )
            {
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

                fvm.source() = rDeltaT *
                    (
                    -coeffs[1] * rho.oldTime().internalField() * vf.oldTime().internalField() * mesh().V0()
                    - coeffs[2] * rho.oldTime().oldTime().internalField() * vf.oldTime().oldTime().internalField() * mesh().V00()
                    - coeffs[3] * rho.oldTime().oldTime().oldTime().internalField() * vf.oldTime().oldTime().oldTime().internalField() * V000.dimensionedInternalField()
                    );
            }
            else
            {
                fvm.source() = rDeltaT * mesh().V() *
                    (
                    -coeffs[1] * rho.oldTime().internalField() * vf.oldTime().internalField()
                    - coeffs[2] * rho.oldTime().oldTime().internalField() * vf.oldTime().oldTime().internalField()
                    - coeffs[3] * rho.oldTime().oldTime().oldTime().internalField() * vf.oldTime().oldTime().oldTime().internalField()
                    );
            }

            return tfvm;
        }

        template<class Type>
        tmp<typename bdf3DdtScheme<Type>::fluxFieldType>
        bdf3DdtScheme<Type>::fvcDdtPhiCorr
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

            scalarField coeffs = coeffs_( U );

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
                    rDeltaT * ddtPhiCoeff
                    * (
                        (
                            (
                                -coeffs[1] * phi.oldTime()
                                - coeffs[2] * phi.oldTime().oldTime()
                                - coeffs[3] * phi.oldTime().oldTime().oldTime()
                            )
                            - (
                                fvc::interpolate
                                (
                                    (
                                        -coeffs[1] * U.oldTime()
                                        - coeffs[2] * U.oldTime().oldTime()
                                        - coeffs[3] * U.oldTime().oldTime().oldTime()
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
        tmp<surfaceScalarField> bdf3DdtScheme<vector>::fvcDdtPhiCorr
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

            scalarField coeffs = coeffs_( U );

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
                const volScalarField & V000 = mesh().objectRegistry::lookupObject<volScalarField>( "V000" );

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

                // Set volume ratio of oldoldoldV over V
                volScalarField V000oV
                (
                    IOobject
                    (
                        "V000oV",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimless,
                    zeroGradientFvPatchScalarField::typeName
                );

                V000oV.internalField() = V000.dimensionedInternalField() / mesh().V();
                V000oV.correctBoundaryConditions();

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
                    surfaceVectorField UV000oV = -fvc::interpolate( V000oV * U.oldTime().oldTime().oldTime() ); // term coming from H/A

                    // Create return variable
                    tmp<surfaceScalarField> ddtPhiCorr
                    (
                        new surfaceScalarField
                        (
                            ddtIOobject,
                            rDeltaT * ddtPhiCoeff
                            * (
                                -coeffs[1] * ( ( UV0oV + fvc::interpolate( V0oV ) * ( Uf.oldTime() ) ) & mesh().Sf() )
                                - coeffs[2] * ( ( UV00oV + fvc::interpolate( V00oV ) * ( Uf.oldTime().oldTime() ) ) & mesh().Sf() )
                                - coeffs[3] * ( ( UV000oV + fvc::interpolate( V000oV ) * ( Uf.oldTime().oldTime().oldTime() ) ) & mesh().Sf() )
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
                    surfaceScalarField phi000 = phi.oldTime().oldTime().oldTime() - ( fvc::interpolate( U.oldTime().oldTime().oldTime() ) & mesh().Sf() );

                    // Create return variable
                    return tmp<fluxFieldType>
                           (
                        new fluxFieldType
                        (
                            ddtIOobject,
                            rDeltaT * ddtPhiCoeff
                            * (
                                -coeffs[1] * phi0
                                - coeffs[2] * phi00
                                - coeffs[3] * phi000
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
        tmp<typename bdf3DdtScheme<Type>::fluxFieldType>
        bdf3DdtScheme<Type>::fvcDdtPhiCorr
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
                    "bdf3DdtScheme<Type>::fvcDdtPhiCorr"
                ) << "dimensions of phiAbs are not correct"
                  << abort( FatalError );

                return fluxFieldType::null();
            }
        }

        template<class Type>
        tmp<surfaceScalarField> bdf3DdtScheme<Type>::meshPhi
            ( const GeometricField<Type, fvPatchField, volMesh> & vf )
        {
            /*scalar deltaT = deltaT_();
             * scalar deltaT0 = deltaT0_(vf);
             *
             * // Coefficient for t-3/2 (between times 0 and 00)
             * scalar coefft0_00 = deltaT/(deltaT + deltaT0);
             *
             * // Coefficient for t-1/2 (between times n and 0)
             * scalar coefftn_0 = 1 + coefft0_00;*/

            scalarField coeffs = coeffs_( vf );
            scalar coefftn_0 = coeffs[0];
            scalar coefft0_00 = coeffs[1] + coeffs[0];
            scalar coefft00_000 = -coeffs[3];

            return coefftn_0 * mesh().phi() + coefft0_00 * mesh().phi().oldTime() + coefft00_000 * mesh().phi().oldTime().oldTime();
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace fv

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
