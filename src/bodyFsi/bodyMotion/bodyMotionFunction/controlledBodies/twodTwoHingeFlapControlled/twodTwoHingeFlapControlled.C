#include "twodTwoHingeFlapControlled.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

using namespace Foam;

defineTypeNameAndDebug(twodTwoHingeFlapControlled, 0);
addToRunTimeSelectionTable(bodyMotionFunction, twodTwoHingeFlapControlled, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > twodTwoHingeFlapControlled::calculatePosition(const scalar time)
{
	tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

	//Determine angle based on force
	scalar rotAngle=0;
	tensor RzCur(cos(rotAngle), -sin(rotAngle), 0, sin(rotAngle), cos(rotAngle), 0, 0, 0, 1);

	//If rotation is postive rotate around top
	vector rotationOrigin = vector::zero;
	if(rotAngle >= 0)
	{
		rotationOrigin = topRotationOrigin_;
	}else{
		rotationOrigin = bottomRotationOrigin_;
	}

	forAll(patchIDs_,ipatch){
		vectorField rotationField = (RzCur & ( initialPoints_[ipatch] - rotationOrigin )) + rotationOrigin;
		position[ipatch] = rotationField;
	}
	return tPosition;
}

void twodTwoHingeFlapControlled::setPatchIDs()
{
	//First fill patchNames_
	patchNames_[0] = patchNameTop_;
	patchNames_[1] = patchNameBottom_;

	forAll(patchNamesBack_,ipatch){
		patchNames_[ipatch+2] = patchNamesBack_[ipatch];
	}

	Info << "patchNames_ = " << patchNames_ << endl;

	forAll(patchNames_,ipatch){
		label patchID = mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]);
		if(patchID < 0){
			FatalErrorIn("void twodTwoHingeFlapControlled::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}

	Info << "patchIDs_ = " << patchIDs_ << endl;
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > twodTwoHingeFlapControlled::getChildMotion()
{
	if(mesh_.time().timeIndex() > localTimeIndex()){
		//Save currentPosition to oldPosition
		forAll(oldPosition_,ipatch){
			oldPosition_[ipatch] = currentPosition_[ipatch];
		}

		//Calculate currentposition based on time
		currentPosition_ = calculatePosition(mesh_.time().value());
	}
	
	tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    forAll(motion,ipatch){
    	motion[ipatch] = currentPosition_[ipatch] - oldPosition_[ipatch];
    }

	return tMotion;
}

tmp<Field<vectorField> > twodTwoHingeFlapControlled::getChildMotion(const Field<vectorField>& patchForces)
{
	if(mesh_.time().timeIndex() > localTimeIndex()){
		//Save currentPosition to oldPosition
		forAll(oldPosition_,ipatch){
			oldPosition_[ipatch] = currentPosition_[ipatch];
		}

		//Calculate currentposition based on time
		currentPosition_ = calculatePosition(mesh_.time().value());
	}
	
	tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    forAll(motion,ipatch){
    	motion[ipatch] = currentPosition_[ipatch] - oldPosition_[ipatch];
    }

	return tMotion;
}

void twodTwoHingeFlapControlled::setInitialPoints(){
	Info << "void twodTwoHingeFlapControlled::setInitialPoints()" << endl;
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList twodTwoHingeFlapControlled::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList twodTwoHingeFlapControlled::getChildMotionPatchNames() const{

	return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twodTwoHingeFlapControlled::twodTwoHingeFlapControlled(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNameTop_(dict.lookup("toppatch")),
patchNameBottom_(dict.lookup("bottompatch")),
patchNamesBack_(dict.lookup("patches")),
patchNames_(patchNamesBack_.size() + 2),
patchIDs_(patchNames_.size(),-1),
proportionalGain_(readScalar(dict.lookup("proportionalGain_"))),
topRotationOrigin_(dict.lookup("topRotationOrigin")),
bottomRotationOrigin_(dict.lookup("bottomRotationOrigin")),
currentPosition_(patchIDs_.size()),
oldPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size())
{
	Info << "twodTwoHingeFlapControlled::twodTwoHingeFlapControlled(const dictionary& dict,const fvMesh& mesh, const word name)" << endl;

	//Set patchIDs_ according to patchNames_
	setPatchIDs();
	
	//Initialize old and current position
	forAll(patchIDs_,ipatch){
		currentPosition_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
		oldPosition_[ipatch] = currentPosition_[ipatch];
	}

	//Set initial points to mesh points
	setInitialPoints();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twodTwoHingeFlapControlled::~twodTwoHingeFlapControlled()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& twodTwoHingeFlapControlled::getPosition() const
{
	return currentPosition_;
}

const Field<vectorField>& twodTwoHingeFlapControlled::getOldPosition() const
{
	return oldPosition_;
}
