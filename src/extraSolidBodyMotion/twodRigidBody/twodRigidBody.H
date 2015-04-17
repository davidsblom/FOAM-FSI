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
    Foam::solidBodyMotionFunctions::twodRigidBody

Description
    Linear oscillation motion function

SourceFiles
    twodRigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef twodRigidBody_H
#define twodRigidBody_H

#include "solidBodyMotionFunction.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                      Class twodRigidBody Declaration
\*---------------------------------------------------------------------------*/

class twodRigidBody
:
    public solidBodyMotionFunction
{
    // Private data

        //Kinematic variables
        scalar rotationAmplitude_;
        scalar rotationFrequency_;
        vector rotationOrigin_;
        scalar rotationPhaseShift_;
        scalar translationAmplitude_;
        scalar translationFrequency_;
        vector translationDirection_;
        scalar translationPhaseShift_;
        bool smoothStart_;

        //Writing properties
        mutable OFstream ofBodyProperties_;
        mutable label ofBodyTimeIndex_;

    // Private Member Functions

        //- Disallow copy construct
        twodRigidBody(const twodRigidBody&);

        //- Disallow default bitwise assignment
        void operator=(const twodRigidBody&);

        //- Write header of properties file
        void writeHeader() const;

        //- Write state of body in properties file
        void writeState() const;

        //- Calculate rotation angle
        scalar calcRotationAngle(const scalar t) const;

        //- Calculate translation vector
        vector calcTranslationVector(const scalar t) const;

        //- Calculate position
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("twodRigidBody");


    // Constructors

        //- Construct from components
        twodRigidBody
        (
            const dictionary& dict,
            const Time& runTime
        );


    // Destructor

        virtual ~twodRigidBody();


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
