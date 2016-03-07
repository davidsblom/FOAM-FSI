#include "airfoilDoubleHingeFlap.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "globalMeshData.H"
#include "processorFvPatch.H"

using namespace Foam;

defineTypeNameAndDebug(airfoilDoubleHingeFlap, 0);
addToRunTimeSelectionTable(bodyMotionFunction, airfoilDoubleHingeFlap, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void airfoilDoubleHingeFlap::writeState(const scalar& time)
{
    scalar flapAngle = getFlapAngle(time)*180.0/mathematicalConstant::pi;
    vector wingState = getWingState(time);//(x,y,theta)
    scalar wingRotAngle = wingState[2];
    vector wingTranslationVector = vector(wingState.x(),wingState.y(),0.0);
    ofBody_ << time << "\t" << wingTranslationVector.x() << "\t" << wingTranslationVector.y() << "\t" << wingRotAngle << "\t" << flapAngle << endl;
}

scalar airfoilDoubleHingeFlap::getFlapAngle(const scalar time) const
{
    const scalar& pi = mathematicalConstant::pi;

    //Determine smoothing value for start
    scalar smoothingScalar = 1;
    if(smoothStart_)
    {
        if(time < 1.0/flapRotationFrequency_){
            smoothingScalar = 0.5 - 0.5*Foam::cos(pi*flapRotationFrequency_*time);

        }
    }

    //Calculate flap angle
    scalar rotAngle = smoothingScalar*flapRotationAmplitude_*sin(2*pi*flapRotationFrequency_*time+flapPhaseShift_);

    return rotAngle;
}

vector airfoilDoubleHingeFlap::getWingState(const scalar& time) const
{
	const scalar& pi = mathematicalConstant::pi;

	scalar smoothingScalarTrans = 1;
	scalar smoothingScalarRot = 1;
	if(smoothStart_)
	{
		if(time < 1.0/translationFrequency_){
			smoothingScalarTrans = 0.5 - 0.5*Foam::cos(pi*translationFrequency_*time);
			//Info << "smoothingScalarTrans = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarTrans*translationAmplitude_ << endl;
		}
		if(time < 1.0/rotationFrequency_){
			smoothingScalarRot = 0.5 - 0.5*Foam::cos(pi*rotationFrequency_*time);
			//Info << "smoothingScalarRot = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarRot*rotationAmplitude_ << endl;
		}
	}

	scalar rotAngle = smoothingScalarRot*rotationAmplitude_*sin(2*pi*rotationFrequency_*time+rotationPhaseShift_);

	vector translationVector =
		smoothingScalarTrans*translationAmplitude_*sin(2*pi*translationFrequency_*time+translationPhaseShift_)*translationDirection_;

	vector state(translationVector.x(),translationVector.y(),rotAngle);

	return state;
}


tmp<Field<vectorField> > airfoilDoubleHingeFlap::calculatePosition(const scalar time)
{
    tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();
    const scalar& pi = mathematicalConstant::pi;

    // ==== First calculate the deformation of the wing (including flap) ==== //
    vector wingState = getWingState(time);//(x,y,theta)
    scalar wingRotAngle = wingState[2];
    vector wingTranslationVector = vector(wingState.x(),wingState.y(),0.0);

    tensor RzCurWing(cos(wingRotAngle), -sin(wingRotAngle), 0, sin(wingRotAngle), cos(wingRotAngle), 0, 0, 0, 1);
    forAll(patchIDs_,ipatch)//do this on all pathces (including flap patches)
    {
        vectorField rotationField = (RzCurWing & ( initialPoints_[ipatch] - rotationOrigin_ )) + rotationOrigin_;
        position[ipatch] = wingTranslationVector + rotationField;
    }

    // === Now calculate additional flap deflection === //
    scalar flapRotAngle = getFlapAngle(time);

    //Take negative angle, since flap angles are defined positive downward (cw)
    tensor RzCurFlap(cos(-flapRotAngle), -sin(-flapRotAngle), 0, sin(-flapRotAngle), cos(-flapRotAngle), 0, 0, 0, 1);

    //Extra transition between rotation point location
    scalar eps_alpha=flapRotationAmplitude_*sin(pi*transitionFraction_);
    scalar Ttrans=transitionFraction_/flapRotationFrequency_;

    //If rotation is postive rotate around top, otherwise around bottom. In the middle do something else
    vector flapRotationOrigin = vector::zero;
    if(flapRotAngle<eps_alpha)
    {
        flapRotationOrigin = topRotationOrigin_;
    }
    else if(flapRotAngle>-eps_alpha)
    {
        flapRotationOrigin = bottomRotationOrigin_;
    }
    else
    {
        scalar trel=asin(flapRotAngle/rotationAmplitude_)/(2.0*pi*rotationFrequency_);//relative t for flap angle
        scalar tnorm=2.0*trel/Ttrans+0.5;//going from 0 to 1

        // flapRotationOrigin=topRotationOrigin_-(topRotationOrigin_-bottomRotationOrigin_)*(tnorm-0.5/pi*(sin(2*pi*tnorm)));
        //x(i) = x0 + dx*(10*localt.^3 - 15*localt.^4 + 6*localt.^5);
        flapRotationOrigin = topRotationOrigin_ - (topRotationOrigin_-bottomRotationOrigin_) * (10*pow(tnorm,3) - 15*pow(tnorm,4) + 6*pow(tnorm,5));
    }

    //Displace rotation point with wing motion
    flapRotationOrigin = (RzCurWing & ( flapRotationOrigin - rotationOrigin_ )) + rotationOrigin_;
    flapRotationOrigin += wingTranslationVector;

    //Save top and bottom position before deforming flap
    vectorField topRotationPoints(topRotationPointIDs_.size(),vector::zero);
    forAll(topRotationPointIDs_,ipoint)
    {
        topRotationPoints[ipoint] = position[0][topRotationPointIDs_[ipoint]];
    }
    vectorField bottomRotationPoints(bottomRotationPointIDs_.size(),vector::zero);
    forAll(bottomRotationPointIDs_,ipoint)
    {
        bottomRotationPoints[ipoint] = position[1][bottomRotationPointIDs_[ipoint]];
    }

    //Set motion to all patches
    forAll(flapPatchIDs_,ipatch)//flapPatches are the first in patchIDs_ and thus in position list
    {
        vectorField rotationField = (RzCurFlap & ( position[ipatch] - flapRotationOrigin )) + flapRotationOrigin;
        position[ipatch] = rotationField;
    }

    //CHECK!Ensure that position of the rotation point remains the same
    forAll(topRotationPointIDs_,ipoint){
        position[0][topRotationPointIDs_[ipoint]] = topRotationPoints[ipoint];
    }
    forAll(bottomRotationPointIDs_,ipoint){
        position[1][bottomRotationPointIDs_[ipoint]] = bottomRotationPoints[ipoint];
    }

    return tPosition;
}

void airfoilDoubleHingeFlap::setFlapRotationPoints()
{
    //Find top rotation point
    Info << "Start searching for top rotation point for flap" << endl;
    topRotationPointIDs_ = findRotationPoints(wingPatchIDs_,flapPatchIDs_[0]);
    const labelList& topFlapPointIDs = mesh_.boundaryMesh()[flapPatchIDs_[0]].meshPoints();
    forAll(topRotationPointIDs_,ipoint){
        if(topRotationPointIDs_[ipoint] >= 0){
            topRotationOrigin_ += mesh_.points()[topFlapPointIDs[topRotationPointIDs_[ipoint]]];
        }else{
            //If values are not set (equals -1), make sure these labelList have size zero. Based on these list the rotation point displacements are set to zero
            //This ensures that this is only done on the processors where these points live.
            topRotationPointIDs_.setSize(0);
        }
    }
    topRotationOrigin_ *= 0.5;
    reduce(topRotationOrigin_,sumOp<vector>());
    Info << "topRotationOrigin = " << topRotationOrigin_ << endl;

    //Find bottom rotation point
    Info << "Start searching for bottom rotation point for flap" << endl;
    bottomRotationPointIDs_ = findRotationPoints(wingPatchIDs_,flapPatchIDs_[1]);
    const labelList& bottomFlapPointIDs = mesh_.boundaryMesh()[flapPatchIDs_[1]].meshPoints();
    forAll(bottomRotationPointIDs_,ipoint){
        if(bottomRotationPointIDs_[ipoint] >= 0){
            bottomRotationOrigin_ += mesh_.points()[bottomFlapPointIDs[bottomRotationPointIDs_[ipoint]]];
        }else{
            //If values are not set (equals -1), make sure these labelList have size zero. Based on these list the rotation point displacements are set to zero
            //This ensures that this is only done on the processors where these points live.
            bottomRotationPointIDs_.setSize(0);
        }
    }
    bottomRotationOrigin_ *= 0.5;
    reduce(bottomRotationOrigin_,sumOp<vector>());
    Info << "bottomRotationOrigin = " << bottomRotationOrigin_ << endl;
}

labelList airfoilDoubleHingeFlap::findRotationPoints(const labelList wingPatchIDs,const label flapPatchID) const
{
    //Prepare label lists
    labelList pointsRotValue(mesh_.nPoints(),0);
    labelList nFoundPoints(Pstream::nProcs(),0);

    const polyPatch& flapPatch(mesh_.boundaryMesh()[flapPatchID]);

    //Loop over wing patches and set value=2 for these points
    forAll(wingPatchIDs,ipatch){
        const label& patchID = wingPatchIDs[ipatch];
        const polyPatch& wingPatch(mesh_.boundaryMesh()[patchID]);
        forAll(wingPatch.meshPoints(),ipoint){
            pointsRotValue[wingPatch.meshPoints()[ipoint]]=2;
        }
    }

    //Loop over flapPatch
    label counter=0;
    labelList rotationPointIDs(2,-1);
    const labelList& patchPointIDs = flapPatch.meshPoints();

       //Add 1 for all flap patch points
    forAll(patchPointIDs,ipoint){
        pointsRotValue[patchPointIDs[ipoint]] += 1;

        //Check if points are marked with 3 (both wing and flap point->rot point)
        if(pointsRotValue[patchPointIDs[ipoint]]==3){
            rotationPointIDs[counter] = ipoint;//Local patch point ID
            counter++;
        }
        if(counter==2){
            break;
        }
    }

    //Check whether the right amount of points are found
    nFoundPoints[Pstream::myProcNo()] = counter;
    reduce(nFoundPoints, sumOp<labelList>());

    //In case the rotation point is located on processor patch
    if(sum(nFoundPoints)!=2){
        globalMeshData gmD(mesh_);

        List<labelList> procPointsValue(gmD.processorPatches().size(),labelList(0,0));
        forAll(gmD.processorPatches(),iprocpatch){
            const processorFvPatch& procPatch(refCast<const processorFvPatch>(mesh_.boundary()[gmD.processorPatches()[iprocpatch]]));

            procPointsValue[iprocpatch] = labelList(mesh_.boundaryMesh()[gmD.processorPatches()[iprocpatch]].meshPoints().size(),0);
            const labelList& procPointIDs = mesh_.boundaryMesh()[gmD.processorPatches()[iprocpatch]].meshPoints();

            //Set value for flap points
            forAll(procPointIDs,ipoint){
                label value(pointsRotValue[procPointIDs[ipoint]]);
                procPointsValue[iprocpatch][ipoint] = value;
            }

            //Send and receive values to/from otherside
            procPatch.send(Pstream::blocking,procPointsValue[iprocpatch]);
            labelList procPointsValueNei(procPointsValue[iprocpatch].size(),0);
            procPatch.receive(Pstream::blocking,procPointsValueNei);

            //Update point values with processorPatch neighbour values
            forAll(procPointsValueNei,ipoint){
                pointsRotValue[procPointIDs[ipoint]] += procPointsValueNei[ipoint];
            }
        }

        //Check now the flap points for value==3 (wing and flap point -> rotation point)
        forAll(flapPatch.meshPoints(),ipoint){
            const label& pointID = flapPatch.meshPoints()[ipoint];
            if(pointsRotValue[pointID]==3){
                rotationPointIDs[counter] = ipoint;//Local patch point ID
                counter++;
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
    }

    return rotationPointIDs;
}

void airfoilDoubleHingeFlap::setFlapPatchIDs()
{
    //Reorder patch names into [top,bottom,rest]
    wordList patchNames(flapPatches_.size());
    patchNames[0] = flapPatches_[patchIndTop_];
    patchNames[1] = flapPatches_[patchIndBottom_];

    forAll(flapPatches_,ipatch){
        if(ipatch != patchIndTop_ && ipatch != patchIndBottom_)
        {
            patchNames[ipatch+2] = flapPatches_[ipatch];
        }
    }
    flapPatches_ = patchNames;

    //Find all patchIDs_
    forAll(flapPatches_,ipatch){
        label patchID = mesh_.boundaryMesh().findPatchID(flapPatches_[ipatch]);
        if(patchID < 0){
            FatalErrorIn("void airfoilDoubleHingeFlap::setPatchIDs()")
            << "patchName " << flapPatches_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            flapPatchIDs_[ipatch] = patchID;
        }
    }
}

void airfoilDoubleHingeFlap::setWingPatchIDs(){
	forAll(wingPatches_,ipatch){
		label patchID(mesh_.boundaryMesh().findPatchID(wingPatches_[ipatch]));
		if(patchID < 0){
			FatalErrorIn("void twodRigidBody::setPatchIDs()")
		    << "patchName " << wingPatches_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			wingPatchIDs_[ipatch] = patchID;
		}
	}
}

void airfoilDoubleHingeFlap::setPatchNamesAndIDs(){
    forAll(flapPatchIDs_,ipatch)
    {
        patchNames_[ipatch] = flapPatches_[ipatch];
        patchIDs_[ipatch] = flapPatchIDs_[ipatch];
    }

    forAll(wingPatchIDs_,ipatch)
    {
        patchNames_[ipatch + flapPatchIDs_.size()] = wingPatches_[ipatch];
        patchIDs_[ipatch + flapPatchIDs_.size()] = wingPatchIDs_[ipatch];
    }
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > airfoilDoubleHingeFlap::getChildMotion()
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

void airfoilDoubleHingeFlap::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
    }
}

const labelList airfoilDoubleHingeFlap::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList airfoilDoubleHingeFlap::getChildMotionPatchNames() const{

    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

airfoilDoubleHingeFlap::airfoilDoubleHingeFlap(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
flapPatches_(dict.lookup("flapPatches")),
flapPatchIDs_(flapPatches_.size()),
patchIndTop_(readLabel(dict.lookup("toppatch"))),
patchIndBottom_(readLabel(dict.lookup("bottompatch"))),
wingPatches_(dict.lookup("wingPatches")),
wingPatchIDs_(wingPatches_.size()),
patchNames_(flapPatchIDs_.size() + wingPatchIDs_.size()),
patchIDs_(flapPatchIDs_.size() + wingPatchIDs_.size()),
rotationAmplitude_(readScalar(dict.lookup("rotationAmplitude"))),
rotationFrequency_(readScalar(dict.lookup("rotationFrequency"))),
rotationOrigin_(dict.lookup("rotationOrigin")),
rotationPhaseShift_(0),
translationAmplitude_(readScalar(dict.lookup("translationAmplitude"))),
translationFrequency_(readScalar(dict.lookup("translationFrequency"))),
translationDirection_(dict.lookup("translationDirection")),
translationPhaseShift_(0),
flapRotationAmplitude_(readScalar(dict.lookup("flapRotationAmplitude"))),
flapRotationFrequency_(readScalar(dict.lookup("flapRotationFrequency"))),
topRotationOrigin_(vector::zero),
bottomRotationOrigin_(vector::zero),
topRotationPointIDs_(2,-1),
bottomRotationPointIDs_(2,-1),
flapPhaseShift_(0),
transitionFraction_(dict.lookupOrDefault("transitionFraction",0.01)),
smoothStart_(false),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    //Wing settings
    translationDirection_ /= mag(translationDirection_) + SMALL;
    rotationAmplitude_ = rotationAmplitude_/180*mathematicalConstant::pi;
    if(dict.found("translationPhaseShift"))
	{
		//Asume shift to be in degrees
		translationPhaseShift_ = readScalar(dict.lookup("translationPhaseShift"))/180.0*mathematicalConstant::pi;
	}
	if(dict.found("rotationPhaseShift"))
	{
		rotationPhaseShift_ = readScalar(dict.lookup("rotationPhaseShift"))/180.0*mathematicalConstant::pi;
	}
	setWingPatchIDs();

    //Flap settings
    flapRotationAmplitude_ = flapRotationAmplitude_/180*mathematicalConstant::pi;
    if(dict.found("phaseShift"))
    {
        flapPhaseShift_ = readScalar(dict.lookup("flapPhaseShift"))/180*mathematicalConstant::pi;
    }
    if(dict.found("smoothStart"))
    {
        smoothStart_ = readBool(dict.lookup("smoothStart"));
    }

    //Set patchIDs_ according to patchNames_
    setFlapPatchIDs();

    //Fill variables with all moving patches
    setPatchNamesAndIDs();

    //Initialize old and current position
    forAll(patchIDs_,ipatch){
        currentPosition_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
        prevPosition_[ipatch] = currentPosition_[ipatch];
    }

    //Set initial points to mesh points
    setInitialPoints();
    setFlapRotationPoints();

    //Prepare body state file if needed. Otherwise remove
    if(writeToFile())
    {
        ofBody_ << "time[s]\tx[m]\ty[m]\ttheta[deg]\tbeta[deg]" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

airfoilDoubleHingeFlap::~airfoilDoubleHingeFlap()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& airfoilDoubleHingeFlap::getPosition() const
{
    return currentPosition_;
}

void airfoilDoubleHingeFlap::update()
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

void airfoilDoubleHingeFlap::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}
