#include "larsFlap.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam;

defineTypeNameAndDebug(larsFlap, 0);
addToRunTimeSelectionTable(bodyMotionFunction, larsFlap, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > larsFlap::calculatePosition(const scalar time)
{
	tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

    //Determine flap angle
	const scalar& pi = mathematicalConstant::pi;
	scalar rotAngle = rotationAmplitude_*sin(2*pi*rotationFrequency_*time+phaseShift_);

	//Constuct rotation matrix to rotate back to normal x,y space
	tensor RzBackward(cos(-reverseAngle_), -sin(-reverseAngle_), 0, sin(-reverseAngle_), cos(-reverseAngle_), 0, 0, 0, 1);

	forAll(patchIDs_,ipatch){

		//Construct displaced point in flap direction space
		vectorField flapDirPoints = initialFlapDirPoints_[ipatch];
		scalarField yflapdir = initialFlapDirPoints_[ipatch].component(1) - sqr(initialFlapDirPoints_[ipatch].component(0))*rotAngle/flapLength_;
		flapDirPoints.replace(1,yflapdir);

		//Rotate back to x,y space
		position[ipatch] = (RzBackward & (flapDirPoints)) + flapCenter_;
	}
	return tPosition;
}

void larsFlap::setPatchIDs()
{
	//First fill patchNames_
	patchNames_[0] = patchNameTop_;
	patchNames_[1] = patchNameBottom_;

	forAll(patchNamesBack_,ipatch){
		patchNames_[ipatch+2] = patchNamesBack_[ipatch];
	}

	forAll(patchNames_,ipatch){
		label patchID = mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]);
		if(patchID < 0){
			FatalErrorIn("void larsFlap::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}
}

void larsFlap::setFlapRotationPoints(){

    //Determine wingPatchIDs
    labelList wingPatchIDs(wingPatches_.size(),-1);
    forAll(wingPatchIDs,ipatch){
    	label patchID = mesh_.boundaryMesh().findPatchID(wingPatches_[ipatch]);
    	if(patchID < 0){
			FatalErrorIn("void larsFlap::setFlapRotationPoints()")
		    << "patchName " << wingPatches_[ipatch] << " in wingPatches does not exist."
		    << abort(FatalError);
		}else{
			wingPatchIDs[ipatch] = patchID;
		}
    }

    //Do top flap patch
    Info << "Start searching for top rotation point for flap" << endl;
    topRotationOrigin_ = findRotationPoints(wingPatchIDs,patchIDs_[0]);
    Info << "topRotationOrigin = " << topRotationOrigin_ << endl;

    Info << "Start searching for bottom rotation point for flap" << endl;
    bottomRotationOrigin_ = findRotationPoints(wingPatchIDs,patchIDs_[1]);
    Info << "bottomRotationOrigin = " << bottomRotationOrigin_ << endl;
}

point larsFlap::findRotationPoints(const labelList wingPatchIDs,const label flapPatchID) const
{
    //Prepare boolean list
    List<bool> pointsBool(mesh_.nPoints(),false);
    labelList nFoundPoints(Pstream::nProcs(),0);

    const polyPatch& flapPatch(mesh_.boundaryMesh()[flapPatchID]);
    forAll(flapPatch.meshPoints(),ipoint){
        pointsBool[flapPatch.meshPoints()[ipoint]]=true;
    }
    //Loop over wingPatchIDs
    label counter=0;
    vector rotationPoint(vector::zero);
    forAll(wingPatchIDs,ipatch){
        const label patchID = wingPatchIDs[ipatch];
        const labelList& patchPointIDs = mesh_.boundaryMesh()[patchID].meshPoints();

        //Check if points are marked
        forAll(patchPointIDs,ipoint){
            if(pointsBool[patchPointIDs[ipoint]]==true){
                rotationPoint += mesh_.points()[patchPointIDs[ipoint]];
                counter++;
            }
            if(counter==2){
                break;
            }
        }
        if(counter==2){
            break;
        }
    }

    //Check whether the right amount of points are found
    nFoundPoints[Pstream::myProcNo()] = counter;
    reduce(nFoundPoints, sumOp<labelList>());
    if(sum(nFoundPoints)!=2){
        FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::setFlapRotationPoints()")
            << "Did not find top rotation point for flap or found to many points. Number of points found = " << sum(nFoundPoints) 
            << abort(FatalError);
    }

    reduce(rotationPoint,sumOp<vector>());
    rotationPoint *=0.5;//Average found points (should be 2 points)
    return rotationPoint;
}

void larsFlap::setFlapDeformationCharacteristics(){

	//Determine flapvector, reverse angle and flapcenter for displacement calculation
	flapVector_ = topRotationOrigin_ - bottomRotationOrigin_;
    flapVector_ /= mag(flapVector_);
    Info << "flapVector_ = " << flapVector_ << endl;
	reverseAngle_ = atan2(flapVector_.x(),flapVector_.y());
	Info << "reverseAngle = " << reverseAngle_ << endl;
	flapCenter_ = 0.5*(topRotationOrigin_+bottomRotationOrigin_);

	//Construct rotation matrices
	tensor RzForward(cos(reverseAngle_), -sin(reverseAngle_), 0, sin(reverseAngle_), cos(reverseAngle_), 0, 0, 0, 1);
	tensor RzBackward(cos(-reverseAngle_), -sin(-reverseAngle_), 0, sin(-reverseAngle_), cos(-reverseAngle_), 0, 0, 0, 1);

	//Determine flap legnth and set initialFlapDirPoints_
	forAll(patchIDs_,ipatch){
		initialFlapDirPoints_[ipatch] = RzForward & (initialPoints_[ipatch] - flapCenter_);
		flapLength_ = max(flapLength_,max(initialFlapDirPoints_[ipatch].component(0)));
	}
	reduce(flapLength_,maxOp<scalar>());//Make sure flapLength is equal to maximum accross processores
	Info << "flapLength = " << flapLength_ << endl;
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > larsFlap::getChildMotion()
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

void larsFlap::setInitialPoints(){
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList larsFlap::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList larsFlap::getChildMotionPatchNames() const{

	return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

larsFlap::larsFlap(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNameTop_(dict.lookup("toppatch")),
patchNameBottom_(dict.lookup("bottompatch")),
patchNamesBack_(dict.lookup("patches")),
patchNames_(patchNamesBack_.size() + 2),
patchIDs_(patchNames_.size(),-1),
wingPatches_(dict.lookup("wingPatches")),
rotationAmplitude_(readScalar(dict.lookup("rotationAmplitude"))),
rotationFrequency_(readScalar(dict.lookup("rotationFrequency"))),
topRotationOrigin_(vector::zero),
bottomRotationOrigin_(vector::zero),
flapVector_(vector::zero),
flapCenter_(vector::zero),
reverseAngle_(0),
flapLength_(0),
phaseShift_(0),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
initialFlapDirPoints_(patchIDs_.size())
{
	rotationAmplitude_ = rotationAmplitude_/180*mathematicalConstant::pi;
	if(dict.found("phaseShift"))
	{
		phaseShift_ = readScalar(dict.lookup("phaseShift"))/180*mathematicalConstant::pi;
	}

	//Set patchIDs_ according to patchNames_
	setPatchIDs();
	
	//Initialize old and current position
	forAll(patchIDs_,ipatch){
		currentPosition_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
		prevPosition_[ipatch] = currentPosition_[ipatch];
	}

	//Set initial points to mesh points
	setInitialPoints();

	//Set flapVector based on rotationPoints
	setFlapRotationPoints();

	//With defined flapVector determine the 
	setFlapDeformationCharacteristics();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

larsFlap::~larsFlap()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& larsFlap::getPosition() const
{
	return currentPosition_;
}

void larsFlap::update()
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
