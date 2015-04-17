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

\*---------------------------------------------------------------------------*/

#include "linearRotation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "tensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearRotation, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearRotation,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidBodyMotionFunctions::linearRotation::writeHeader() const
{
    ofBodyProperties_ << "time [s]\tangle [deg]" << endl;
}

void Foam::solidBodyMotionFunctions::linearRotation::writeState() const
{
    scalar t = time_.value();

    if(time_.timeIndex()>ofBodyTimeIndex_){
        const scalar rotAngle = calcRotationAngle(t)*180/mathematicalConstant::pi;
        ofBodyProperties_ << t << "\t" << rotAngle << endl;

        ofBodyTimeIndex_ = time_.timeIndex();
    }
}

Foam::scalar Foam::solidBodyMotionFunctions::linearRotation::calcRotationAngle(const scalar t) const
{
    const scalar& pi = mathematicalConstant::pi;
    scalar rotAngle = 2*pi*rotationFrequency_*t;

    return rotAngle;
}

Foam::septernion
Foam::solidBodyMotionFunctions::linearRotation::calcTransformation(const scalar t) const
{
    const scalar rotAngle = calcRotationAngle(t);

    quaternion R(0,0,rotAngle);
    septernion TR(septernion(rotationOrigin_)*R*septernion(-rotationOrigin_));

    return TR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearRotation::linearRotation
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    rotationFrequency_(readScalar(SBMFCoeffs_.lookup("rotationFrequency"))),
    rotationOrigin_(SBMFCoeffs_.lookup("rotationOrigin")),
    ofBodyProperties_(type()+"-properties.dat"),
    ofBodyTimeIndex_(runTime.timeIndex())
{
    writeHeader();
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearRotation::~linearRotation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::linearRotation::transformation() const
{
    scalar t = time_.value();

    septernion TR = calcTransformation(t);

    writeState();

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::linearRotation::velocity() const
{
    // Velocity is calculated as a difference
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    return (calcTransformation(t + dt) - calcTransformation(t))/dt;
}


bool Foam::solidBodyMotionFunctions::linearRotation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("rotationFrequency") >> rotationFrequency_;
    SBMFCoeffs_.lookup("rotationOrigin") >> rotationOrigin_;
    
    return true;
}


// ************************************************************************* //
