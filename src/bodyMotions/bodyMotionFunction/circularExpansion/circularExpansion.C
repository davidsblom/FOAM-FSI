#include "circularExpansion.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

using namespace Foam;

defineTypeNameAndDebug(circularExpansion, 0);
addToRunTimeSelectionTable(bodyMotionFunction, circularExpansion, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void circularExpansion::writeState(const scalar& time)
{
    scalar state = getState(time);
    ofBody_ << time << "\t" << state << endl;
}

scalar circularExpansion::getState(const scalar& time)
{
    const scalar& pi = mathematicalConstant::pi;

    scalar radius = radius0_ + amplitude_*0.5*(1.0-cos(2*pi*frequency_*time));

    return radius;
}

tmp<Field<vectorField> > circularExpansion::calculatePosition(const scalar time)
{

    tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

    scalar state = getState(time);

    forAll(patchIDs_,ipatch){
        position[ipatch] = initialPoints_[ipatch];
        vectorField r = (initialPoints_[ipatch] - axis_*(axis_ & initialPoints_[ipatch]));
        r -= origin_;
        forAll(r,ipoint)
        {
            scalar radius = mag(r[ipoint]);
            if(radius<state){
                position[ipatch][ipoint] += (state-radius)*(r[ipoint]/radius);
            }
        }
    }


    return tPosition;
}

void circularExpansion::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void circularExpansion::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > circularExpansion::getChildMotion()
{
    //If the localTimeIndex is not equal to mesh time index, update has not been called
    if(mesh_.time().timeIndex() != localTimeIndex())
    {
        WarningIn("circularExpansion::getChildMotion()")
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

void circularExpansion::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.Cf().boundaryField()[patchIDs_[ipatch]];
    }
}

const labelList circularExpansion::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList circularExpansion::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

circularExpansion::circularExpansion(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
radius0_(readScalar(dict.lookup("initialRadius"))),
amplitude_(readScalar(dict.lookup("amplitude"))),
frequency_(readScalar(dict.lookup("frequency"))),
origin_(dict.lookup("origin")),
axis_(dict.lookup("axis")),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    Info << "A = " << amplitude_ << ", origin = " << origin_ << endl;

    if(mag(axis_)>SMALL)
    {
        axis_/=mag(axis_);
    }
    else
    {
        FatalErrorIn("void circularExpansion::Constructor()")
        << "mag(axis) == 0. Give axis with magnitude. Obviously ;)"
        << abort(FatalError);
    }

    //Set patchIDs_ according to patchNames_
    setPatchIDs();

    //Initialize old and current position
    forAll(patchIDs_,ipatch){
        currentPosition_[ipatch] = mesh_.Cf().boundaryField()[patchIDs_[ipatch]];
        prevPosition_[ipatch] = currentPosition_[ipatch];

    }

    //Set initial points to mesh points
    setInitialPoints();

    if(writeToFile())
    {
        ofBody_ << "time [s]\tr [m]" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

circularExpansion::~circularExpansion()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& circularExpansion::getPosition() const
{
    return currentPosition_;
}

void circularExpansion::update()
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

void circularExpansion::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}
