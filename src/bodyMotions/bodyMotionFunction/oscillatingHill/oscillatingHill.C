#include "oscillatingHill.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "foamTime.H"

using namespace Foam;

defineTypeNameAndDebug(oscillatingHill, 0);
addToRunTimeSelectionTable(bodyMotionFunction, oscillatingHill, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void oscillatingHill::writeState(const scalar& time)
{
    scalar state = getState(time);
    ofBody_ << time << "\t" << state << endl;
}

scalar oscillatingHill::getState(const scalar& time)
{
    const scalar& pi = mathematicalConstant::pi;

    scalar state = amplitude_*0.5*(1.0-cos(2*pi*frequency_*time));

    return state;
}

tmp<Field<vectorField> > oscillatingHill::calculatePosition(const scalar time)
{
    tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

    const scalar& pi = mathematicalConstant::pi;
    scalar state = getState(time);

    forAll(patchIDs_,ipatch)
    {
        position[ipatch] = initialPoints_[ipatch];
        forAll(position[ipatch],ipoint)
        {
            scalar Xt1=((initialPoints_[ipatch][ipoint]-origin_) & tangent1_);
            scalar Yt2=((initialPoints_[ipatch][ipoint]-origin_) & tangent2_);

            position[ipatch][ipoint] += state*direction_*(0.5-0.5*cos(2*pi*(Xt1-0.5)))*(0.5-0.5*cos(2*pi*(Yt2-0.5)));
        }
    }


    return tPosition;
}

void oscillatingHill::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void oscillatingHill::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

void oscillatingHill::setTangents()
{
    tangent1_=(vector(direction_.y(),direction_.z(),direction_.x()) ^ direction_);
    tangent2_=(direction_ ^ tangent1_);

    Info << "tangent1_ = " << tangent1_ << endl;
    Info << "tangent2_ = " << tangent2_ << endl;

    /*Info << "Check normal and tangents:"
    << " n x t1 = " << (direction_ & tangent1_)
    << " n x t2 = " << (direction_ & tangent2_)
    << " t1 x t2 = " << (tangent1_ & tangent2_) << endl;*/
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > oscillatingHill::getChildMotion()
{
    //If the localTimeIndex is not equal to mesh time index, update has not been called
    if(mesh_.time().timeIndex() != localTimeIndex())
    {
        WarningIn("oscillatingHill::getChildMotion()")
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

void oscillatingHill::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.Cf().boundaryField()[patchIDs_[ipatch]];
    }
}

const labelList oscillatingHill::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList oscillatingHill::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillatingHill::oscillatingHill(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
amplitude_(readScalar(dict.lookup("amplitude"))),
frequency_(readScalar(dict.lookup("frequency"))),
origin_(dict.lookup("origin")),
direction_(dict.lookup("direction")),
tangent1_(vector::zero),
tangent2_(vector::zero),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    direction_ /=(mag(direction_)+SMALL);
    setTangents();
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

oscillatingHill::~oscillatingHill()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& oscillatingHill::getPosition() const
{
    return currentPosition_;
}

void oscillatingHill::update()
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

void oscillatingHill::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}
