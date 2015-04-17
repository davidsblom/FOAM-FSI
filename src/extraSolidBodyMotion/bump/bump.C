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

#include "bump.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(bump, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        bump,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::solidBodyMotionFunctions::bump::calcPosition
(
    const scalar t
) const
{
    /*scalar currentScale = 0.0;
    if(t<=period_){
        currentScale =  0.5*(1.0 - cos(2*pi*t/period_));
    }
    
    return amplitude_*currentScale;*/

    scalar currentScale = 1.0;
    if(t<=period_){
        currentScale = 0.5*(t - period_/pi/2*sin(2*pi/period_*t));
    }else{
        currentScale = 0.5*period_;
    }
    
    return amplitude_*currentScale;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::bump::bump
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    amplitude_(SBMFCoeffs_.lookup("amplitude")),
    period_(readScalar(SBMFCoeffs_.lookup("period")))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::bump::~bump()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::bump::transformation() const
{
    scalar t = time_.value();

    septernion TR(calcPosition(t), quaternion::I);

    Info<< "solidBodyMotionFunctions::bump::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::bump::velocity() const
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


bool Foam::solidBodyMotionFunctions::bump::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("amplitude") >> amplitude_;
    SBMFCoeffs_.lookup("amplitude") >> period_;

    return true;
}


// ************************************************************************* //
