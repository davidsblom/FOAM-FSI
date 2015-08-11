#include "twodRigidBody.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "foamTime.H"

using namespace Foam;

defineTypeNameAndDebug(twodRigidBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, twodRigidBody, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void twodRigidBody::writeState(const scalar& time)
{
    vector state = getState(time);
    ofBody_ << time << "\t" << state.z()*180.0/mathematicalConstant::pi << "\t" << state.x() << "\t" << state.y()<< endl;
}

vector twodRigidBody::getState(const scalar& time)
{
    const scalar& pi = mathematicalConstant::pi;

    scalar smoothingScalarTrans = 1;
    scalar smoothingScalarRot = 1;
    if(smoothStart_)
    {
        if(time < 1.0/translationFrequency_){
            smoothingScalarTrans = 0.5 - 0.5*Foam::cos(pi*translationFrequency_*time);
            Info << "smoothingScalarTrans = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarTrans*translationAmplitude_ << endl;
        }
        if(time < 1.0/rotationFrequency_){
            smoothingScalarRot = 0.5 - 0.5*Foam::cos(pi*rotationFrequency_*time);
            Info << "smoothingScalarRot = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarRot*rotationAmplitude_ << endl;
        }
    }

    scalar rotAngle = smoothingScalarRot*rotationAmplitude_*sin(2*pi*rotationFrequency_*time+translationPhaseShift_+phaseShift_);

    vector translationVector =
        smoothingScalarTrans*translationAmplitude_*sin(2*pi*translationFrequency_*time+translationPhaseShift_)*translationDirection_;

    vector state(translationVector.x(),translationVector.y(),rotAngle);

    return state;
}

tmp<Field<vectorField> > twodRigidBody::calculatePosition(const scalar time)
{
    tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

    vector state = getState(time);
    scalar rotAngle = state.z();
    vector translationVector(state.x(),state.y(),0.0);

    /*const scalar& pi = mathematicalConstant::pi;

    scalar rotAngle = rotationAmplitude_*sin(2*pi*rotationFrequency_*time+translationPhaseShift_+phaseShift_);

    vector translationVector =
        translationAmplitude_*sin(2*pi*translationFrequency_*time+translationPhaseShift_)*translationDirection_;*/

    tensor RzCur(cos(rotAngle), -sin(rotAngle), 0, sin(rotAngle), cos(rotAngle), 0, 0, 0, 1);
    forAll(patchIDs_,ipatch){
        vectorField rotationField = (RzCur & ( initialPoints_[ipatch] - rotationOrigin_ )) + rotationOrigin_;
        position[ipatch] = translationVector + rotationField;
    }

    return tPosition;
}

void twodRigidBody::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void twodRigidBody::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > twodRigidBody::getChildMotion()
{
    //If the localTimeIndex is not equal to mesh time index, update has not been called
    if(mesh_.time().timeIndex() != localTimeIndex())
    {
        WarningIn("twodRigidBody::getChildMotion()")
            << "Local time index is not equal to global time index. Update functions needs to be called" << endl;
        update();
    }

    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    forAll(motion,ipatch){
        //Set motion
        motion[ipatch] = currentPosition_[ipatch] - prevPosition_[ipatch];
        //Update prevPosition
        prevPosition_[ipatch] = currentPosition_[ipatch];
    }

    return tMotion;
}

void twodRigidBody::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        //initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
        initialPoints_[ipatch] = mesh_.Cf().boundaryField()[patchIDs_[ipatch]];
    }
}

const labelList twodRigidBody::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList twodRigidBody::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twodRigidBody::twodRigidBody(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
rotationAmplitude_(readScalar(dict.lookup("rotationAmplitude"))),
rotationFrequency_(readScalar(dict.lookup("rotationFrequency"))),
translationAmplitude_(readScalar(dict.lookup("translationAmplitude"))),
translationFrequency_(readScalar(dict.lookup("translationFrequency"))),
translationDirection_(dict.lookup("translationDirection")),
translationPhaseShift_(0),
rotationOrigin_(dict.lookup("rotationOrigin")),
phaseShift_(0),
smoothStart_(false),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{

    translationDirection_ /= mag(translationDirection_) + SMALL;
    rotationAmplitude_ = rotationAmplitude_/180.0*mathematicalConstant::pi;
    if(dict.found("translationPhaseShift"))
    {
        //Asume shift to be in degrees
        translationPhaseShift_ = readScalar(dict.lookup("translationPhaseShift"))/180.0*mathematicalConstant::pi;
    }
    if(dict.found("phaseShift"))
    {
        phaseShift_ = readScalar(dict.lookup("phaseShift"))/180.0*mathematicalConstant::pi;
    }
    if(dict.found("smoothStart"))
    {
        smoothStart_ = readBool(dict.lookup("smoothStart"));
    }

    //Set patchIDs_ according to patchNames_
    setPatchIDs();

    //Initialize old and current position
    forAll(patchIDs_,ipatch){
        //currentPosition_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
        currentPosition_[ipatch] = mesh_.Cf().boundaryField()[patchIDs_[ipatch]];
        prevPosition_[ipatch] = currentPosition_[ipatch];

    }

    //Set initial points to mesh points
    setInitialPoints();

    if(writeToFile())
    {
        ofBody_ << "time [s]\ttheta [deg]\tx [m]\ty [m]" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twodRigidBody::~twodRigidBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& twodRigidBody::getPosition() const
{
    return currentPosition_;
}

void twodRigidBody::update()
{
    if(localTimeIndex() < mesh_.time().timeIndex()){
        forAll(prevPosition_,ipatch){
            prevPosition_[ipatch] = currentPosition_[ipatch];
        }
    }
    bodyMotionFunction::update();//This updates localTimeIndex

    //Since this is a forced motion the position can be calculated here
    currentPosition_ = calculatePosition(mesh_.time().value());
}

void twodRigidBody::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}
