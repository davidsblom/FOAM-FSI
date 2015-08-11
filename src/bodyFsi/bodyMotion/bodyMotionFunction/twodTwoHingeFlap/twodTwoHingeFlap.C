#include "twodTwoHingeFlap.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "globalMeshData.H"
#include "processorFvPatch.H"

using namespace Foam;

defineTypeNameAndDebug(twodTwoHingeFlap, 0);
addToRunTimeSelectionTable(bodyMotionFunction, twodTwoHingeFlap, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void twodTwoHingeFlap::writeState(const scalar& time)
{
    scalar flapAngle = getFlapAngle(time)*180.0/mathematicalConstant::pi;
    ofBody_ << time << "\t" << flapAngle << endl;
}

scalar twodTwoHingeFlap::getFlapAngle(const scalar time)
{
    const scalar& pi = mathematicalConstant::pi;

    //Determine smoothing value for start
    scalar smoothingScalar = 1;
    if(smoothStart_)
    {
        if(time < 1.0/rotationFrequency_){
            smoothingScalar = 0.5 - 0.5*Foam::cos(pi*rotationFrequency_*time);
            Info << "smoothingScalar = " << smoothingScalar << ", amplitude = " << smoothingScalar*rotationAmplitude_ << endl;

        }
    }

    //Calculate flap angle
    scalar rotAngle = smoothingScalar*rotationAmplitude_*sin(2*pi*rotationFrequency_*time+phaseShift_);

    return rotAngle;
}

tmp<Field<vectorField> > twodTwoHingeFlap::calculatePosition(const scalar time)
{
    tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

    const scalar& pi = mathematicalConstant::pi;
    scalar rotAngle = getFlapAngle(time);

    //Take negative angle, since flap angles are defined positive downward (cw)
    tensor RzCur(cos(-rotAngle), -sin(-rotAngle), 0, sin(-rotAngle), cos(-rotAngle), 0, 0, 0, 1);

    //Extra transition between rotation point location
    scalar eps=0.01;
    scalar eps_alpha=rotationAmplitude_*sin(2*pi*eps);
    scalar Ttrans=2.0*eps*1.0/rotationFrequency_;

    //If rotation is postive rotate around top
    vector rotationOrigin = vector::zero;
    if(rotAngle<-eps_alpha)
    {
        rotationOrigin = topRotationOrigin_;
    }
    else if(rotAngle>eps_alpha)
    {
        rotationOrigin = bottomRotationOrigin_;
    }
    else
    {
        scalar trel=asin(rotAngle/rotationAmplitude_)/(2.0*pi*rotationFrequency_);//relative t for flap angle
        //scalar tnorm=trel/Ttrans;//going from -0.5 to 0.5
        //rotationOrigin=topRotationOrigin_-(topRotationOrigin_-bottomRotationOrigin_)*(0.5+0.5*(sin(pi*tnorm)));
        scalar tnorm=trel/Ttrans+0.5;//going from 0 to 1
        rotationOrigin=topRotationOrigin_-(topRotationOrigin_-bottomRotationOrigin_)*(tnorm-0.5/pi*(sin(2*pi*tnorm)));
    }
/*	if(rotAngle <= 0)
    {
        rotationOrigin = topRotationOrigin_;
    }else{
        rotationOrigin = bottomRotationOrigin_;
    }*/

    //Set motion to all patches
    forAll(patchIDs_,ipatch){
        vectorField rotationField = (RzCur & ( initialPoints_[ipatch] - rotationOrigin )) + rotationOrigin;
        position[ipatch] = rotationField;
    }

    //Ensure that position of the rotation point remains the same
    forAll(topRotationPointIDs_,ipoint){
        position[0][topRotationPointIDs_[ipoint]] = initialPoints_[0][topRotationPointIDs_[ipoint]];
    }
    forAll(bottomRotationPointIDs_,ipoint){
        position[1][bottomRotationPointIDs_[ipoint]] = initialPoints_[1][bottomRotationPointIDs_[ipoint]];
    }

    return tPosition;
}

void twodTwoHingeFlap::setPatchIDs()
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
            FatalErrorIn("void twodTwoHingeFlap::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

void twodTwoHingeFlap::setFlapRotationPoints(){

    //Determine wingPatchIDs
    labelList wingPatchIDs(wingPatches_.size(),-1);
    forAll(wingPatchIDs,ipatch){
        label patchID = mesh_.boundaryMesh().findPatchID(wingPatches_[ipatch]);
        if(patchID < 0){
            FatalErrorIn("void twodTwoHingeFlap::setFlapRotationPoints()")
            << "patchName " << wingPatches_[ipatch] << " in wingPatches does not exist."
            << abort(FatalError);
        }else{
            wingPatchIDs[ipatch] = patchID;
        }
    }

    //Find top rotation point
    Info << "Start searching for top rotation point for flap" << endl;
    topRotationPointIDs_ = findRotationPoints(wingPatchIDs,patchIDs_[0]);
    const labelList& topFlapPointIDs = mesh_.boundaryMesh()[patchIDs_[0]].meshPoints();
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
    bottomRotationPointIDs_ = findRotationPoints(wingPatchIDs,patchIDs_[1]);
    const labelList& bottomFlapPointIDs = mesh_.boundaryMesh()[patchIDs_[1]].meshPoints();
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

labelList twodTwoHingeFlap::findRotationPoints(const labelList wingPatchIDs,const label flapPatchID) const
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

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > twodTwoHingeFlap::getChildMotion()
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

void twodTwoHingeFlap::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
    }
}

const labelList twodTwoHingeFlap::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList twodTwoHingeFlap::getChildMotionPatchNames() const{

    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twodTwoHingeFlap::twodTwoHingeFlap(const dictionary& dict,const fvMesh& mesh, const word name):
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
topRotationPointIDs_(2,-1),
bottomRotationPointIDs_(2,-1),
phaseShift_(0),
smoothStart_(false),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    rotationAmplitude_ = rotationAmplitude_/180*mathematicalConstant::pi;
    if(dict.found("phaseShift"))
    {
        phaseShift_ = readScalar(dict.lookup("phaseShift"))/180*mathematicalConstant::pi;
    }
    if(dict.found("smoothStart"))
    {
        smoothStart_ = readBool(dict.lookup("smoothStart"));
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
    setFlapRotationPoints();

    //Prepare body state file if needed. Otherwise remove
    if(writeToFile())
    {
        ofBody_ << "time[s]\tbeta[deg]" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twodTwoHingeFlap::~twodTwoHingeFlap()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& twodTwoHingeFlap::getPosition() const
{
    return currentPosition_;
}

void twodTwoHingeFlap::update()
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

void twodTwoHingeFlap::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}
