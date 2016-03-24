
#ifndef bdf2DdtScheme_H
#define bdf2DdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace fv
    {
        /*---------------------------------------------------------------------------*\
        *                      Class bdf2DdtScheme Declaration
        \*---------------------------------------------------------------------------*/

        template<class Type>
        class bdf2DdtScheme
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

            // - Disallow default bitwise copy construct
            bdf2DdtScheme( const bdf2DdtScheme & );

            // - Disallow default bitwise assignment
            void operator=( const bdf2DdtScheme & );

            public:
                // - Runtime type information
                TypeName( "bdf2" );


                // Constructors

                // - Construct from mesh
                bdf2DdtScheme( const fvMesh & mesh )
                    :
                    ddtScheme<Type>( mesh )
                {}

                // - Construct from mesh and Istream
                bdf2DdtScheme(
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
        tmp<surfaceScalarField> bdf2DdtScheme<scalar>::fvcDdtPhiCorr(
            const volScalarField & rA,
            const volScalarField & U,
            const surfaceScalarField & phi
            );


        template<>
        tmp<surfaceScalarField> bdf2DdtScheme<scalar>::fvcDdtPhiCorr(
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
 #   include "bdf2DdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
