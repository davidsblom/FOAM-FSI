
#ifndef bdf3DdtScheme_H
#define bdf3DdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace fv
    {
        /*---------------------------------------------------------------------------*\
        *                      Class bdf3DdtScheme Declaration
        \*---------------------------------------------------------------------------*/

        template<class Type>
        class bdf3DdtScheme
            :
            public fv::ddtScheme<Type>
        {
            // Private Member Functions

            // - Return the current time-step
            scalar deltaT_() const;

            // - Return the previous time-step
            scalar deltaT0_() const;

            // - Return the previous time-step or GREAT if the old timestep field
            // wasn't available in which case Euler ddt is used
            template<class GeoField>
            scalar deltaT0_( const GeoField & ) const;

            // - Return the coefficients for integration. Needed for first time steps
            template<class GeoField>
            scalarField coeffs_( const GeoField & ) const;

            // - Return the standard coefficients for integration.
            scalarField coeffs_() const;

            // - Disallow default bitwise copy construct
            bdf3DdtScheme( const bdf3DdtScheme & );

            // - Disallow default bitwise assignment
            void operator=( const bdf3DdtScheme & );

public:

            // - Runtime type information
            TypeName( "bdf3" );


            // Constructors

            // - Construct from mesh
            bdf3DdtScheme( const fvMesh & mesh )
                :
                ddtScheme<Type>( mesh )
            {}

            // - Construct from mesh and Istream
            bdf3DdtScheme(
                const fvMesh & mesh,
                Istream & is
                )
                :
                ddtScheme<Type>( mesh, is )
            {}


            // Member Functions

            // - Return mesh reference
            const fvMesh & mesh() const
            {
                return fv::ddtScheme<Type>::mesh();
            }

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt( const dimensioned<Type> & );

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt( const GeometricField<Type, fvPatchField, volMesh> & );

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt(
                const dimensionedScalar &,
                const GeometricField<Type, fvPatchField, volMesh> &
                );

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt(
                const volScalarField &,
                const GeometricField<Type, fvPatchField, volMesh> &
                );

            tmp<fvMatrix<Type> > fvmDdt( GeometricField<Type, fvPatchField, volMesh> & );

            tmp<fvMatrix<Type> > fvmDdt(
                const dimensionedScalar &,
                GeometricField<Type, fvPatchField, volMesh> &
                );

            tmp<fvMatrix<Type> > fvmDdt(
                const volScalarField &,
                GeometricField<Type, fvPatchField, volMesh> &
                );

            typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

            tmp<fluxFieldType> fvcDdtPhiCorr(
                const volScalarField & rA,
                const GeometricField<Type, fvPatchField, volMesh> & U,
                const fluxFieldType & phi
                );

            tmp<fluxFieldType> fvcDdtPhiCorr(
                const volScalarField & rA,
                const volScalarField & rho,
                const GeometricField<Type, fvPatchField, volMesh> & U,
                const fluxFieldType & phi
                );


            tmp<surfaceScalarField> meshPhi( const GeometricField<Type, fvPatchField, volMesh> & );
        };


        template<>
        tmp<surfaceScalarField> bdf3DdtScheme<scalar>::fvcDdtPhiCorr(
            const volScalarField & rA,
            const volScalarField & U,
            const surfaceScalarField & phi
            );


        template<>
        tmp<surfaceScalarField> bdf3DdtScheme<scalar>::fvcDdtPhiCorr(
            const volScalarField & rA,
            const volScalarField & rho,
            const volScalarField & U,
            const surfaceScalarField & phi
            );


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace fv

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
 #   include "bdf3DdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
