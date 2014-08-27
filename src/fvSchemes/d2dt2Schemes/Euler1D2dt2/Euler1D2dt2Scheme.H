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

Class
    Foam::fv::Euler1D2dt2Scheme

Description
    First-order Euler implicit d2dt2 using the current and two previous
    time-step values.

SourceFiles
    Euler1D2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef Euler1D2dt2Scheme_H
#define Euler1D2dt2Scheme_H

#include "d2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class Euler1D2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Euler1D2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        Euler1D2dt2Scheme(const Euler1D2dt2Scheme&);

        //- Disallow default bitwise assignment
        void operator=(const Euler1D2dt2Scheme&);


public:

    //- Runtime type information
    TypeName("Euler1");


    // Constructors

        //- Construct from mesh
        Euler1D2dt2Scheme(const fvMesh& mesh)
        :
            d2dt2Scheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        Euler1D2dt2Scheme(const fvMesh& mesh, Istream& is)
        :
            d2dt2Scheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Euler1D2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
