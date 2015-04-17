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

#include "smoothTranslation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(smoothTranslation, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        smoothTranslation,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::solidBodyMotionFunctions::smoothTranslation::calcPosition
(
    const scalar t
) const
{
    scalar currentScale = 1.0;
    if(t<=dT_){
        currentScale = 0.5*(t - dT_/pi*sin(pi/dT_*t));
    }else{
        currentScale = t - 0.5*dT_;
    }
    
    return amplitude_*currentScale;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::smoothTranslation::smoothTranslation
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    amplitude_(SBMFCoeffs_.lookup("amplitude")),
    dT_(readScalar(SBMFCoeffs_.lookup("dT")))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::smoothTranslation::~smoothTranslation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::smoothTranslation::transformation() const
{
    scalar t = time_.value();

    septernion TR(calcPosition(t), quaternion::I);

    Info<< "solidBodyMotionFunctions::smoothTranslation::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::smoothTranslation::velocity() const
{
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    septernion TV
    (
        (calcPosition(t + dt) - calcPosition(t))/dt,
        quaternion::zero
    );

    return TV;
}


bool Foam::solidBodyMotionFunctions::smoothTranslation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("amplitude") >> amplitude_;
    SBMFCoeffs_.lookup("dT") >> dT_;

    return true;
}


// ************************************************************************* //
