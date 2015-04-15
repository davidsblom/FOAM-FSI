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
    Foam::SRF::rpmAbs

Description
    Basic SRF model whereby angular velocity is specified in terms of
    a (global) axis and revolutions-per-minute [rpmAbs]

SourceFiles
    rpmAbs.C

\*---------------------------------------------------------------------------*/

#ifndef SRFModelAbsrpmAbs_H
#define SRFModelAbsrpmAbs_H

#include "SRFModelAbs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SRF
{

/*---------------------------------------------------------------------------*\
                            Class rpmAbs Declaration
\*---------------------------------------------------------------------------*/

class rpmAbs
:
    public SRFModelAbs
{

    // Private data

        //- Revolutions per minute
        scalar rpm_;


    // Private member functions

        //- Disallow default bitwise copy construct
        rpmAbs(const rpmAbs&);

        //- Disallow default bitwise assignment
        void operator=(const rpmAbs&);


public:

    //- Runtime type information
    TypeName("rpm");


    // Constructors

        //- Construct from components
        rpmAbs(const volVectorField& U);


    // Destructor

        ~rpmAbs();

    // Member functions

        // I-O

            //- Read
            bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SRF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
