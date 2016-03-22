
#ifndef bdf1DdtScheme_H
#define bdf1DdtScheme_H

#include "ddtScheme.H"
#include "volFields.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace fv
    {
        /*---------------------------------------------------------------------------*\
        *                      Class bdf1DdtScheme Declaration
        \*---------------------------------------------------------------------------*/

        template<class Type>
        class bdf1DdtScheme
              :
                public ddtScheme<Type>
        {
            // Private Member Functions

            // - Disallow default bitwise copy construct
            bdf1DdtScheme( const bdf1DdtScheme & );

            // - Disallow default bitwise assignment
            void operator=( const bdf1DdtScheme & );

            public:
                // - Runtime type information
                TypeName( "bdf1" );


                // Constructors

                // - Construct from mesh
                bdf1DdtScheme( const fvMesh & mesh )
                    :
                    ddtScheme<Type>( mesh )
                {}

                // - Construct from mesh and Istream
                bdf1DdtScheme(
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
        tmp<surfaceScalarField> bdf1DdtScheme<scalar>::fvcDdtPhiCorr(
            const volScalarField & rA,
            const volScalarField & U,
            const surfaceScalarField & phi
            );


        template<>
        tmp<surfaceScalarField> bdf1DdtScheme<scalar>::fvcDdtPhiCorr(
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
 #   include "bdf1DdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
