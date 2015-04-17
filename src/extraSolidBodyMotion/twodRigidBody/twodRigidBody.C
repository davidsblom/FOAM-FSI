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

#include "twodRigidBody.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "tensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(twodRigidBody, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        twodRigidBody,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidBodyMotionFunctions::twodRigidBody::writeHeader() const
{
    ofBodyProperties_ << "time [s]\tangle [deg]\tx [m]\ty [m]" << endl;
}

void Foam::solidBodyMotionFunctions::twodRigidBody::writeState() const
{
    scalar t = time_.value();

    if(time_.timeIndex()>ofBodyTimeIndex_){
        const scalar rotAngle = calcRotationAngle(t)*180/mathematicalConstant::pi;
        const vector translationVector = calcTranslationVector(t);
        ofBodyProperties_ << t << "\t" << rotAngle << "\t" << translationVector.x() << "\t" << translationVector.y() << endl;

        ofBodyTimeIndex_ = time_.timeIndex();
    }
}

Foam::scalar Foam::solidBodyMotionFunctions::twodRigidBody::calcRotationAngle(const scalar t) const
{
    const scalar& pi = mathematicalConstant::pi;

    //Smoothing if requested
    scalar smoothingScalarRot = 1;
    if(smoothStart_)
    {
        if(t < 1.0/(rotationFrequency_+SMALL)){
            smoothingScalarRot = 0.5 - 0.5*Foam::cos(pi*rotationFrequency_*t);
            //Info << "smoothingScalarRot = " << smoothingScalarRot << ", amplitude = " << smoothingScalarRot*rotationAmplitude_ << endl;
        }
    }

    //Calculate angle
    scalar rotAngle = smoothingScalarRot*rotationAmplitude_*sin(2*pi*rotationFrequency_*t+rotationPhaseShift_);

    return rotAngle;
}

Foam::vector Foam::solidBodyMotionFunctions::twodRigidBody::calcTranslationVector(const scalar t) const
{
    Info << "Foam::vector Foam::solidBodyMotionFunctions::twodRigidBody::calcTranslationVector(const scalar t) const" << endl;
    const scalar& pi = mathematicalConstant::pi;

    //Smoothing if requested
    scalar smoothingScalarTrans = 1;
    if(smoothStart_)
    {
        if(t < 1.0/(translationFrequency_+SMALL)){
            smoothingScalarTrans = 0.5 - 0.5*Foam::cos(pi*translationFrequency_*t);
            //Info << "smoothingScalarTrans = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarTrans*translationAmplitude_ << endl;
        }
    }

    //Calculate vector
    vector translationVector = smoothingScalarTrans*translationAmplitude_*sin(2*pi*translationFrequency_*t+translationPhaseShift_)*translationDirection_;

    return translationVector;
}

Foam::septernion
Foam::solidBodyMotionFunctions::twodRigidBody::calcTransformation(const scalar t) const
{
    scalar rotAngle = calcRotationAngle(t);
    vector translationVector = calcTranslationVector(t);

    quaternion R(0,0,rotAngle);
    septernion TR(septernion(rotationOrigin_ + translationVector)*R*septernion(-rotationOrigin_));

    return TR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::twodRigidBody::twodRigidBody
(
    const dictionary& dict,
    const Time& runTime
)
:
    solidBodyMotionFunction(dict, runTime),
    rotationAmplitude_(readScalar(SBMFCoeffs_.lookup("rotationAmplitude"))),
    rotationFrequency_(readScalar(SBMFCoeffs_.lookup("rotationFrequency"))),
    rotationOrigin_(SBMFCoeffs_.lookup("rotationOrigin")),
    rotationPhaseShift_(0),
    translationAmplitude_(readScalar(SBMFCoeffs_.lookup("translationAmplitude"))),
    translationFrequency_(readScalar(SBMFCoeffs_.lookup("translationFrequency"))),
    translationDirection_(SBMFCoeffs_.lookup("translationDirection")),
    translationPhaseShift_(0),
    smoothStart_(false),
    ofBodyProperties_(type()+"-properties.dat"),
    ofBodyTimeIndex_(runTime.timeIndex())
{
    translationDirection_ /= mag(translationDirection_) + SMALL;
    rotationAmplitude_ = rotationAmplitude_/180*mathematicalConstant::pi;
    Info << "Foam::solidBodyMotionFunctions::twodRigidBody::twodRigidBody(dict,time)" << endl;
    if(SBMFCoeffs_.found("translationPhaseShift"))
    {
        //Asume shift to be in degrees
        translationPhaseShift_ = readScalar(SBMFCoeffs_.lookup("translationPhaseShift"))/180.0*mathematicalConstant::pi;
    }
    if(SBMFCoeffs_.found("rotationPhaseShift"))
    {
        rotationPhaseShift_ = readScalar(SBMFCoeffs_.lookup("rotationPhaseShift"))/180.0*mathematicalConstant::pi;
    }
    if(SBMFCoeffs_.found("smoothStart"))
    {
        smoothStart_ = readBool(SBMFCoeffs_.lookup("smoothStart"));
		Info << "Doing smooth start for first period of motions" << nl << endl;
    }

    writeHeader();
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::twodRigidBody::~twodRigidBody()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::twodRigidBody::transformation() const
{
    scalar t = time_.value();

    septernion TR = calcTransformation(t);

    /*Info<< "solidBodyMotionFunctions::twodRigidBody::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;*/
    writeState();

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::twodRigidBody::velocity() const
{
    // Velocity is calculated as a difference
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    return (calcTransformation(t + dt) - calcTransformation(t))/dt;
}


bool Foam::solidBodyMotionFunctions::twodRigidBody::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("rotationAmplitude") >> rotationAmplitude_;
    SBMFCoeffs_.lookup("rotationFrequency") >> rotationFrequency_;
    SBMFCoeffs_.lookup("rotationOrigin") >> rotationOrigin_;
    SBMFCoeffs_.lookup("translationAmplitude") >> translationAmplitude_;
    SBMFCoeffs_.lookup("translationFrequency") >> translationFrequency_;
    SBMFCoeffs_.lookup("translationDirection") >> translationDirection_;

    return true;
}


// ************************************************************************* //
