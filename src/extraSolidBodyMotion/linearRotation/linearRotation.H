/*---------------------------------------------------------------------------*\
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

Class
    Foam::solidBodyMotionFunctions::linearRotation

Description
    Linear oscillation motion function

SourceFiles
    linearRotation.C

\*---------------------------------------------------------------------------*/

#ifndef linearRotation_H
#define linearRotation_H

#include "solidBodyMotionFunction.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                      Class linearRotation Declaration
\*---------------------------------------------------------------------------*/

class linearRotation
:
    public solidBodyMotionFunction
{
    // Private data

        //Kinematic variables
        scalar rotationFrequency_;
        vector rotationOrigin_;

        //Writing properties
        mutable OFstream ofBodyProperties_;
        mutable label ofBodyTimeIndex_;

    // Private Member Functions

        //- Disallow copy construct
        linearRotation(const linearRotation&);

        //- Disallow default bitwise assignment
        void operator=(const linearRotation&);

        //- Write header of properties file
        void writeHeader() const;

        //- Write state of body in properties file
        void writeState() const;

        //- Calculate rotation angle
        scalar calcRotationAngle(const scalar t) const;

        //- Calculate position
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("linearRotation");


    // Constructors

        //- Construct from components
        linearRotation
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Destructor

        virtual ~linearRotation();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Return the solid-body motion velocity
        virtual septernion velocity() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
