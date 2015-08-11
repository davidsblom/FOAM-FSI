#include "beamForced.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam;

defineTypeNameAndDebug(beamForced, 0);
addToRunTimeSelectionTable(bodyMotionFunction, beamForced, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void beamForced::writeState(const scalar time)
{
	scalar tipDisp = getTipDisplacement(time);
	ofBody_ << time << "\t" <<  tipDisp << endl;
}

scalar beamForced::getTipDisplacement(const scalar time)
{
	const scalar& pi = mathematicalConstant::pi;

	scalar currentscale = 0.0;
	scalar dT_=1.0/frequency_;
	if(time>T0_)
	{
		currentscale = sin(2*pi*frequency_*(time-T0_));
		if(time<T0_+dT_){
			scalar smoothFactor=0.5*(1-cos(pi/dT_*(time-T0_)));
			currentscale *= smoothFactor;
		}
	}

	//Adding something weird with dt dependency for testing
	//currentscale *= (1+10*sqr(mesh_.time().deltaT().value()));

	//Info << "currentScale = " << currentscale << endl;
	return currentscale*amplitude_;
}


tmp<Field<vectorField> > beamForced::calculatePosition(const scalar time)
{
	tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

	scalar tipDisp = getTipDisplacement(time);

	forAll(patchIDs_,ipatch){
		vectorField deltaPoints(initialPoints_.size(),vector::zero);
		deltaPoints = transDir_ * tipDisp*Foam::sqr(((initialPoints_[ipatch] - origin_) & normalDir_)/beamLength_);

		position[ipatch] = initialPoints_[ipatch] + deltaPoints;

		//scalarField x = initialPoints_[ipatch].component(0) - origin_.x();
		//position[ipatch].replace(1, initialPoints_[ipatch].component(1) + 0.28065306122*sqr(x)*currentscale);
	}

	return tPosition;
}

void beamForced::setPatchIDs(){
	forAll(patchNames_,ipatch){
		label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
		if(patchID < 0){
			FatalErrorIn("void beamForced::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > beamForced::getChildMotion()
{
	//If the localTimeIndex is not equal to mesh time index, update has not been called
	if(mesh_.time().timeIndex() != localTimeIndex())
	{
		WarningIn("beamForced::getChildMotion()")
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

void beamForced::setBeamLength(){

	scalar maxPos = -1e15;
	scalar minPos = 1e15;
	forAll(initialPoints_,ipatch){
		scalarField normalInitialPoints = initialPoints_[ipatch] & normalDir_;

		scalar maxLocal = max(normalInitialPoints);
		maxPos = max(maxPos,maxLocal);

		scalar minLocal = min(normalInitialPoints);
		minPos = min(minPos,minLocal);
	}

	beamLength_ = maxPos - minPos;
}

void beamForced::setInitialPoints(){
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList beamForced::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList beamForced::getChildMotionPatchNames() const{
	return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

beamForced::beamForced(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
amplitude_(readScalar(dict.lookup("amplitude"))),
frequency_(readScalar(dict.lookup("frequency"))),
origin_(dict.lookup("origin")),
normalDir_(dict.lookup("normalDirection")),
transDir_(dict.lookup("translationDirection")),
beamLength_(1.0),
T0_(mesh.time().startTime().value()),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
	T0_ = dict.lookupOrDefault("startTime",T0_);
	Info << "Start time of structure = " << T0_ << endl;
	//Set patchIDs_ according to patchNames_
	setPatchIDs();

	//Initialize old and current position
	forAll(patchIDs_,ipatch){
		currentPosition_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
		prevPosition_[ipatch] = currentPosition_[ipatch];
	}

	//Set initial points to mesh points
	setInitialPoints();

	//Calculate beam length
	normalDir_ /= mag(normalDir_);
	transDir_ /= mag(transDir_);
	setBeamLength();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

beamForced::~beamForced()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& beamForced::getPosition() const
{
	return currentPosition_;
}

void beamForced::update()
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

void beamForced::write()
{
	if(writeToFile())
	{
		writeState(mesh_.time().value());
	}
}
