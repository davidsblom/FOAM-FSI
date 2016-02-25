#include "beamRotationForced.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam;

defineTypeNameAndDebug(beamRotationForced, 0);
addToRunTimeSelectionTable(bodyMotionFunction, beamRotationForced, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void beamRotationForced::writeState(const scalar time)
{
	scalar rotAngle = getRotationAngle(time);
	ofBody_ << time << "\t" <<  rotAngle << endl;
}

scalar beamRotationForced::getRotationAngle(const scalar time) const
{
    const scalar& pi = mathematicalConstant::pi;

    scalar rotationAmplitude_ = amplitude_/180.0*pi;
    //Calculate flap angle
    scalar rotAngle = rotationAmplitude_*sin(2*pi*frequency_*time);

    return rotAngle;
}

/*scalar beamRotationForced::getTipDisplacement(const scalar time)
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

	return currentscale*amplitude_;
}*/


/*tmp<Field<vectorField> > beamRotationForced::calculatePosition(const scalar time)
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
}*/

tmp<Field<vectorField> > beamRotationForced::calculatePosition(const scalar time)
{
    tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

    // ==== First calculate the deformation of the wing (including flap) ==== //
    scalar beamRotAngle = getRotationAngle(time);


    /*forAll(patchIDs_,ipatch)//do this on all pathces
    {
        vectorField rotationField = ((initialPoints_[ipatch] - origin_ ) * cos(beamRotAngle)
                                    + (rotationAxis_ ^ ((initialPoints_[ipatch] - origin_ ))) * sin(beamRotAngle)
                                    + rotationAxis_ * (rotationAxis_ & ( initialPoints_[ipatch] - origin_ )) * (1-cos(beamRotAngle)) ) + origin_;
        position[ipatch] = rotationField;
    }

    //Add the bending to the beam
    forAll(patchIDs_,ipatch)//do this on all pathces
    {
        scalarField transDisp = (initialPoints_[ipatch] & transDir_) + ((position[ipatch] - initialPoints_[ipatch]) & transDir_) * Foam::sqr(((initialPoints_[ipatch] - origin_) & normalDir_)/beamLength_);
        scalarField normDisp = ((position[ipatch] - origin_) & normalDir_);
        position[ipatch] = normDisp*normalDir_ + transDisp*transDir_ + origin_;
    }*/

    forAll(patchIDs_,ipatch)//do this on all pathces
    {
        //Get local rotation angle based on theta*(xn/L)^2
        scalarField beamLocalRotAngle = beamRotAngle*Foam::sqr(((initialPoints_[ipatch] - origin_) & normalDir_)/beamLength_);
        vectorField rotationField = ((initialPoints_[ipatch] - origin_ ) * cos(beamLocalRotAngle)
                                    + (rotationAxis_ ^ ((initialPoints_[ipatch] - origin_ ))) * sin(beamLocalRotAngle)
                                    + rotationAxis_ * (rotationAxis_ & ( initialPoints_[ipatch] - origin_ )) * (1.0 - cos(beamLocalRotAngle)) ) + origin_;
        position[ipatch] = rotationField;
    }

    return tPosition;
}

void beamRotationForced::setPatchIDs(){
	forAll(patchNames_,ipatch){
		label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
		if(patchID < 0){
			FatalErrorIn("void beamRotationForced::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > beamRotationForced::getChildMotion()
{
	//If the localTimeIndex is not equal to mesh time index, update has not been called
	if(mesh_.time().timeIndex() != localTimeIndex())
	{
		WarningIn("beamRotationForced::getChildMotion()")
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

void beamRotationForced::setBeamLength(){

	scalar maxPos = -1e15;
	scalar minPos = 1e15;
	forAll(initialPoints_,ipatch){
		scalarField normalInitialPoints = initialPoints_[ipatch] & normalDir_;

		scalar maxLocal = max(normalInitialPoints);
		maxPos = max(maxPos,maxLocal);

		scalar minLocal = min(normalInitialPoints);
		minPos = min(minPos,minLocal);
	}
    //Ensure for parallel calculations the correct length is achieved
    reduce( maxPos, maxOp<scalar>() );
    reduce( minPos, minOp<scalar>() );

	beamLength_ = maxPos - minPos;
}

void beamRotationForced::setInitialPoints(){
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList beamRotationForced::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList beamRotationForced::getChildMotionPatchNames() const{
	return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

beamRotationForced::beamRotationForced(const dictionary& dict,const fvMesh& mesh, const word name):
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
rotationAxis_(dict.lookup("rotationAxis")),
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
    rotationAxis_ /= mag(rotationAxis_);
	setBeamLength();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

beamRotationForced::~beamRotationForced()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& beamRotationForced::getPosition() const
{
	return currentPosition_;
}

void beamRotationForced::update()
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

void beamRotationForced::write()
{
	if(writeToFile())
	{
		writeState(mesh_.time().value());
	}
}
