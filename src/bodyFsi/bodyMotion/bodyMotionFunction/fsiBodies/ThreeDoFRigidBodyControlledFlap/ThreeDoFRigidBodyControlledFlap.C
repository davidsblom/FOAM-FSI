
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    ThreeDoFRigidBodyControlledFlap

\*---------------------------------------------------------------------------*/
#include "ThreeDoFRigidBodyControlledFlap.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"
#include "mathematicalConstants.H"
#include "twoDPointCorrector.H"

using namespace Foam;

defineTypeNameAndDebug(ThreeDoFRigidBodyControlledFlap, 0);
addToRunTimeSelectionTable(bodyMotionFunction, ThreeDoFRigidBodyControlledFlap, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
void ThreeDoFRigidBodyControlledFlap::writeState(const scalar& time)
{
	//Body
    ofBody_ << time << "\t" << bodyState_[0] << "\t" << bodyState_[1] << "\t" << bodyState_[2] 
            << "\t" << bodyState_[3] << "\t" << bodyState_[4] << "\t" << bodyState_[5] 
            << "\t" << bodyForces_[0] << "\t" << bodyForces_[1] << "\t" << bodyForces_[2]
            << "\t" << beta_ << endl;
}

vector ThreeDoFRigidBodyControlledFlap::calculateForcesAndMoments(const Field<vectorField>& patchForces) const
{ 
    vector forcesAndMoments(0,0,0);
    scalar xbody = bodyState_[0];
    scalar ybody = bodyState_[1];
    vector CofR = rc_ + vector(xbody,ybody,0);

    vector totalMoments(vector::zero);
    vector totalForces(vector::zero);
    forAll(patchIDs_,ipatch)
    {
        vectorField Md = mesh_.boundaryMesh()[patchIDs_[ipatch]].faceCentres() - CofR;
        vectorField moments = (Md ^ patchForces[patchIDs_[ipatch]]);

        //Parallel implementation
        totalMoments += gSum(moments);
        totalForces += gSum(patchForces[patchIDs_[ipatch]]);
    }

    Info << "totalForces = " << totalForces << endl;
    Info << "totalMoments = " << totalMoments << endl;

    forcesAndMoments.x()=totalForces.x();
    forcesAndMoments.y()=totalForces.y();
    forcesAndMoments.z()=totalMoments.z();

    return forcesAndMoments;
}

vector ThreeDoFRigidBodyControlledFlap::calculateWingPosition(const Field<vectorField>& patchForces)
{
    const scalar& dt = mesh_.time().deltaT().value();

    vector forcesAndMoments = calculateForcesAndMoments(patchForces);
    bodyForces_ = forcesAndMoments;

    //Start by saving previous iteration
    bodyStatePrevIter_ = bodyState_;

    //Set values of x,y,theta,dxdt,dydt,dthetadt to old values initially (every subiteration of fsi also)
    scalar xbody = bodyStateOld_[0];
    scalar ybody = bodyStateOld_[1];
    scalar theta = bodyStateOld_[2];
    scalar dxbody = bodyStateOld_[3];
    scalar dybody = bodyStateOld_[4];
    scalar dtheta = bodyStateOld_[5];

    if(integrationScheme_=="explicit")
    {
        Info << "Using explicit body motion calculation " << endl;

        simpleMatrix<scalar> C(3,0,0);
        C[0][0] = m_;
        C[0][2] = -m_*Lcg2rc_*sin(theta+thetaGeom_);
        C[1][1] = m_;
        C[1][2] = m_*Lcg2rc_*cos(theta+thetaGeom_);
        C[2][0] = -m_*Lcg2rc_*sin(theta+thetaGeom_);
        C[2][1] = m_*Lcg2rc_*cos(theta+thetaGeom_);
        C[2][2] = Icg_ + m_*pow(Lcg2rc_,2);

        C.source()[0] = forcesAndMoments.x() - c_.x()*dxbody - k_.x()*xbody + m_*Lcg2rc_*pow(dtheta,2)*cos(theta+thetaGeom_);
        C.source()[1] = forcesAndMoments.y() - c_.y()*dybody - k_.y()*ybody + m_*Lcg2rc_*pow(dtheta,2)*sin(theta+thetaGeom_);
        C.source()[2] = forcesAndMoments.z() - ctheta_*dtheta - ktheta_*theta;
        List<scalar> bodyStateAcc = C.LUsolve();

        //scalarList bodyStateVel(6,0.0);
        //bodyStateVel[0]=dxbody;bodyStateVel[1]=dybody;bodyStateVel[2]=dtheta;
        //bodyStateVel[3]=bodyStateAcc[0];bodyStateVel[4]=bodyStateAcc[1];bodyStateVel[5]=bodyStateAcc[2];
        bodyStateVel_[0]=dxbody;bodyStateVel_[1]=dybody;bodyStateVel_[2]=dtheta;
        bodyStateVel_[3]=bodyStateAcc[0];bodyStateVel_[4]=bodyStateAcc[1];bodyStateVel_[5]=bodyStateAcc[2];

        bodyState_ = bodyStateOld_ + dt*bodyStateVel_;
    }
    else if(integrationScheme_=="RK4")
    {
        List<scalar> rkCoeffs(4,1.0);
        rkCoeffs[0]=1.0/4.0;rkCoeffs[1]=1.0/3.0;rkCoeffs[2]=1.0/2.0;
        
        simpleMatrix<scalar> C(3,0,0);
        for(int i=0;i<rkCoeffs.size();i++){
            
            C = simpleMatrix<scalar>(3,0,0);
            C[0][0] = m_;
            C[0][2] = -m_*Lcg2rc_*sin(theta+thetaGeom_);
            C[1][1] = m_;
            C[1][2] = m_*Lcg2rc_*cos(theta+thetaGeom_);
            C[2][0] = -m_*Lcg2rc_*sin(theta+thetaGeom_);
            C[2][1] = m_*Lcg2rc_*cos(theta+thetaGeom_);
            C[2][2] = Icg_ + m_*pow(Lcg2rc_,2);

            C.source()[0] = forcesAndMoments.x() - c_.x()*dxbody - k_.x()*xbody + m_*Lcg2rc_*pow(dtheta,2)*cos(theta+thetaGeom_);
            C.source()[1] = forcesAndMoments.y() - c_.y()*dybody - k_.y()*ybody + m_*Lcg2rc_*pow(dtheta,2)*sin(theta+thetaGeom_);
            C.source()[2] = forcesAndMoments.z() - ctheta_*dtheta - ktheta_*theta;//Add some moment here
            List<scalar> bodyStateAcc = C.LUsolve();

            //scalarList bodyStateVel(6,0.0);
            //bodyStateVel[0]=dxbody;bodyStateVel[1]=dybody;bodyStateVel[2]=dtheta;
            //bodyStateVel[3]=bodyStateAcc[0];bodyStateVel[4]=bodyStateAcc[1];bodyStateVel[5]=bodyStateAcc[2];
        	bodyStateVel_[0]=dxbody;bodyStateVel_[1]=dybody;bodyStateVel_[2]=dtheta;
        	bodyStateVel_[3]=bodyStateAcc[0];bodyStateVel_[4]=bodyStateAcc[1];bodyStateVel_[5]=bodyStateAcc[2];            

            List<scalar> bodyStageState = bodyStateOld_ + rkCoeffs[i]*dt*bodyStateVel_;
            xbody = bodyStageState[0];ybody = bodyStageState[1];theta = bodyStageState[2];
            dxbody = bodyStageState[3];dybody = bodyStageState[4];dtheta = bodyStageState[5];
        }
        
        bodyState_[0]=xbody;
        bodyState_[1]=ybody;
        bodyState_[2]=theta;
        bodyState_[3]=dxbody;
        bodyState_[4]=dybody;
        bodyState_[5]=dtheta;
    }
    else
    {
        FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::calculateWingPosition()")
            << "integrationScheme " << integrationScheme_ << " not available."
            << "Available schemes are: (explicit, RK4)"
            << abort(FatalError);
    }
    
    vector bodyStatePos(bodyState_[0],bodyState_[1],bodyState_[2]);
    Info << "bodyState = " << bodyStatePos << endl;
    
    return bodyStatePos;
}

tmp<Field<vectorField> > ThreeDoFRigidBodyControlledFlap::calculateFlapMotion()
{
    //Return variable
    tmp<Field<vectorField> > tMotion(new Field<vectorField>(flapPatches_.size()));
    Field<vectorField>& motion = tMotion();

    if(flapUpdated_){
        Info << "flap already updated. Beta = " << beta_ << endl;
        forAll(flapPatches_,ipatch){
            motion[ipatch] = 0*initialPoints_[flapPatches_[ipatch]];
        }
        return tMotion;
    }else
    {
        const scalar& time = mesh_.time().value();
        const scalar& dt = mesh_.time().deltaT().value();
        const scalar& pi = mathematicalConstant::pi;

        //Set beta based on y velocity of body
        scalar betaMax = 7.0/180.0*pi;
        if(time>controllerStartTime_){

        	//Controller
			scalar controllerOutput = controllerVelocityGain_*bodyStateOld_[4] + controllerAccGain_*bodyStateVelOld_[4];

        	if(flapIntegrationScheme_=="Euler"){
            	//First order time integration
            	beta_ = betaOld_ + controllerOutput*dt/180.0*pi;    
            	if(mag(beta_)>betaMax){
                	beta_ = sign(beta_)*betaMax;
            	}
            }else if(flapIntegrationScheme_=="BDF2"){
            	//Second order time integration
    			scalar deltaT = mesh_.time().deltaT().value();
    			scalar deltaT0 = mesh_.time().deltaT0().value();
    			if(mesh_.time().timeIndex()==1){
    				deltaT0=GREAT;
    			}

    			scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    			scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    			scalar coefft0  = coefft + coefft00;
				
				beta_ = coefft0/coefft*betaOld_ - coefft00/coefft*betaOldOld_ + 1.0/coefft*(controllerOutput*dt)/180.0*pi;
            	//beta_ = 4.0/3.0*betaOld_ -1.0/3.0*betaOldOld_ + 2.0/3.0*((controllerVelocityGain_)*bodyState_[4]*dt)/180.0*pi;
            	if(mag(beta_)>betaMax){
                	beta_ = sign(beta_)*betaMax;
            	}
            }
        }else{
            beta_ = betaOld_;
        }
        
        //Rotation matrices
        tensor RzCur
        (
            Foam::cos(-beta_), -Foam::sin(-beta_), 0,
            Foam::sin(-beta_),  Foam::cos(-beta_), 0,
            0, 0, 1
        );

        tensor RzOld
        (
            Foam::cos(-betaOld_), -Foam::sin(-betaOld_), 0,
            Foam::sin(-betaOld_),  Foam::cos(-betaOld_), 0,
            0, 0, 1
        );

        Info << "Calculated rotation matrices" << endl;
        //If rotation is postive rotate around top
        vector rotationOrigin = vector::zero;
        if(beta_ >= 0)
        {
			rotationOrigin = bottomRotationOrigin_;
        }else{
            
            rotationOrigin = topRotationOrigin_;
        }

        //If rotation is postive rotate around top
        vector rotationOriginOld = vector::zero;
        if(betaOld_ >= 0){
        	rotationOriginOld = bottomRotationOrigin_;
        }else{
            rotationOriginOld = topRotationOrigin_;
        }
        Info << "Set rotation origins" << endl;

        forAll(flapPatches_,ipatch){
            vectorField rotationField = (RzCur & ( initialPoints_[flapPatches_[ipatch]] - rotationOrigin )) + rotationOrigin;
            vectorField rotationFieldOld = (RzOld & ( initialPoints_[flapPatches_[ipatch]] - rotationOriginOld )) + rotationOriginOld;
            motion[ipatch] = rotationField - rotationFieldOld;
        }

        Info << "Flap position is now updated. Beta = " << beta_ << endl;
        flapUpdated_ = true;

        return tMotion;
    }
}

void ThreeDoFRigidBodyControlledFlap::readInitialStateIfPresent()
{
    vector state = dict_.lookupOrDefault("initialState",vector::zero);
    vector stateVel = dict_.lookupOrDefault("initialStateVelocity",vector::zero);

    bodyState_[0] = state.x();
    bodyState_[1] = state.y();
    bodyState_[2] = state.z();
    bodyState_[3] = stateVel.x();
    bodyState_[4] = stateVel.y();
    bodyState_[5] = stateVel.z();
    theta0_ = bodyState_[2]; //needed in rotation tensors

    bodyStateOldOld_ = bodyStateOld_;
    bodyStateOld_ = bodyState_;

    bodyState0_ = bodyState_;

    Info << "Initial state is set to: " << bodyState_ << endl;
}

void ThreeDoFRigidBodyControlledFlap::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }

    //Check whether flapPatches are within length of patches
    forAll(flapPatches_,ipatch){
        if(flapPatches_[ipatch] > patchIDs_.size())
        {
            FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::setPatchIDs()")
            << "entry " << ipatch << " of flapPatches: " << flapPatches_[ipatch] << " is outside of provided patch list"
            << abort(FatalError);
        }
    }

    //Check whether flapPatchTop is within patches
    if(flapPatchTop_ > patchIDs_.size())
    {
        FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::setPatchIDs()")
        << "flapPatchTop is set to " << flapPatchTop_ << ", which is outside of provided patch list"
        << abort(FatalError);
    }
    if(flapPatchBottom_ > patchIDs_.size())
    {
        FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::setPatchIDs()")
        << "flapPatchBottom is set to " << flapPatchBottom_ << ", which is outside of provided patch list"
        << abort(FatalError);
    }
}

void ThreeDoFRigidBodyControlledFlap::setFlapRotationPoints(){

    //Determine wingPatchIDs
    List<bool> isFlapPatch(patchIDs_.size(),false);
    forAll(flapPatches_,ipatch){
        isFlapPatch[flapPatches_[ipatch]] = true;
    }

    labelList wingPatchIDs(patchIDs_.size(),-1);
    label counter=0;
    forAll(isFlapPatch,ipatch){
        if(isFlapPatch[ipatch]==false){
            wingPatchIDs[counter] = patchIDs_[ipatch];
            counter++;
        }
    }
    wingPatchIDs.setSize(counter);

    //Do top flap patch
    Info << "Start searching for top rotation point for flap" << endl;
    topRotationOrigin_ = findRotationPoints(wingPatchIDs,flapPatchTop_);
    Info << "topRotationOrigin = " << topRotationOrigin_ << endl;

    Info << "Start searching for bottom rotation point for flap" << endl;
    bottomRotationOrigin_ = findRotationPoints(wingPatchIDs,flapPatchBottom_);
    Info << "bottomRotationOrigin = " << bottomRotationOrigin_ << endl;
}

point ThreeDoFRigidBodyControlledFlap::findRotationPoints(const labelList wingPatchIDs,const label flapPatchIndex) const
{
    //Prepare boolean list
    List<bool> pointsBool(mesh_.nPoints(),false);
    labelList nFoundPoints(Pstream::nProcs(),0);

    const polyPatch& flapPatch(mesh_.boundaryMesh()[patchIDs_[flapPatchIndex]]);
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

    //Set correct rotation point
    reduce(rotationPoint,sumOp<vector>());
    rotationPoint *=0.5;//Average found points (should be 2 points)
    return rotationPoint;
}

// * * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

tmp<Field<vectorField> > ThreeDoFRigidBodyControlledFlap::getChildMotion()
{
    FatalErrorIn("void ThreeDoFRigidBodyControlledFlap::getChildMotion()")
            << "this bodyMotionFunction cannot be used in a non-fsi solver, since getChildMotion(const Field<vectorField>& patchForces) needs to be called"
            << abort(FatalError);

    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    //Field<vectorField>& motion = tMotion();

    return tMotion;
}

tmp<Field<vectorField> > ThreeDoFRigidBodyControlledFlap::getChildMotion(const Field<vectorField>& patchForces)
{
    Info << "Solve Solid " << name_ << endl;
    //Prepare return variable
    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    //Calculate the bodyState based on the forces
    calculateWingPosition(patchForces);

    //Calculate position based on bodystate
    scalar thetaCur   = bodyState_[2];
    scalar thetaPrev   = bodyStatePrevIter_[2];
    //It is assumed that rotation center (in the mesh) is located at rc_ + initial translation
    vector CofR = rc_ + vector(bodyState0_[0],bodyState0_[1],0);

    //Translation is subtracting current state with previous state
    vector translationVector(bodyState_[0] - bodyStatePrevIter_[0],bodyState_[1] - bodyStatePrevIter_[1],0);

    //Rotation
    tensor RzCur
    (
        Foam::cos(thetaCur-theta0_), -Foam::sin(thetaCur-theta0_), 0,
        Foam::sin(thetaCur-theta0_),  Foam::cos(thetaCur-theta0_), 0,
        0, 0, 1
    );

    tensor RzPrev
    (
        Foam::cos(thetaPrev-theta0_), -Foam::sin(thetaPrev-theta0_), 0,
        Foam::sin(thetaPrev-theta0_),  Foam::cos(thetaPrev-theta0_), 0,
        0, 0, 1
    );

    forAll(patchIDs_,ipatch){
        vectorField rotationField((RzCur-RzPrev) & (initialPoints_[ipatch] - CofR));
        motion[ipatch] = translationVector + (rotationField);        
    }

    //Set flap position
    Field<vectorField> flapMotion = calculateFlapMotion();
    forAll(flapMotion,ipatch){

        motion[flapPatches_[ipatch]] += flapMotion[ipatch];
    }

    return tMotion;
}

void ThreeDoFRigidBodyControlledFlap::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
    }
}

const labelList ThreeDoFRigidBodyControlledFlap::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList ThreeDoFRigidBodyControlledFlap::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ThreeDoFRigidBodyControlledFlap::ThreeDoFRigidBodyControlledFlap(const dictionary& dict,const fvMesh& mesh,const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
flapPatches_(dict.lookup("flapPatches")),
flapPatchTop_(readLabel(dict.lookup("flapPatchTop"))),
flapPatchBottom_(readLabel(dict.lookup("flapPatchBottom"))),
integrationScheme_(dict_.lookup("integrationScheme")),
flapIntegrationScheme_(dict.lookupOrDefault("flapIntegrationScheme",word("BDF2"))),
k_(dict_.lookup("k")),
c_(dict_.lookup("c")),
ktheta_(readScalar(dict_.lookup("ktheta"))),
ctheta_(readScalar(dict_.lookup("ctheta"))),
m_(0.0),
Icg_(readScalar(dict_.lookup("Icg"))),
Lcg2rc_(readScalar(dict_.lookup("Lcg2rc"))),
rc_(dict_.lookup("rc")),
thetaGeom_(readScalar(dict_.lookup("thetaGeom"))),
topRotationOrigin_(vector::zero),
bottomRotationOrigin_(vector::zero),
controllerStartTime_(readScalar(dict.lookup("controllerStartTime"))),
controllerVelocityGain_(readScalar(dict.lookup("controllerVelocityGain"))),
controllerAccGain_(readScalar(dict.lookup("controllerAccGain"))),
bodyState_(6,0.0),
bodyStatePrevIter_(6,0.0),
bodyStateOld_(6,0.0),
bodyStateOldOld_(6,0.0),
bodyState0_(6,0.0),
bodyStateVel_(6,0.0),
bodyStateVelOld_(6,0.0),
theta0_(0),
beta_(0),
betaOld_(0),
betaOldOld_(0),
flapUpdated_(false),
bodyForces_(vector::zero),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    setPatchIDs();
    setInitialPoints();
    setFlapRotationPoints();

	//Determine thickness of mesh and set m accordingly
	Vector<label> e = (-mesh.geometricD()+Vector<label>(1,1,1))/2;
	scalar meshDepth(vector(e.x(),e.y(),e.z()) & vector(mesh.bounds().span()));
    scalar m(readScalar(dict_.lookup("m")));
	m_ = m*meshDepth;

    //Read initial state
    readInitialStateIfPresent();

    Info << "Using flapIntegrationScheme: " << flapIntegrationScheme_ << endl;

    if(writeToFile())
    {
        ofBody_ << "time\tx\ty\ttheta\tdxdt\tdydt\tdthetadt\tFx\tFy\tMtheta\tbeta" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ThreeDoFRigidBodyControlledFlap::~ThreeDoFRigidBodyControlledFlap()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*vector ThreeDoFRigidBodyControlledFlap::calculateSteadyStatePosition(const vectorField& forces)
{
    vector forcesAndMoments = calculateForcesAndMoments(forces);
    bodyForces_ = forcesAndMoments;

    //Start by saving previous iteration
    bodyStatePrevIter_ = bodyState_;

    //Calculate steady state deformation
    bodyState_[0] = forcesAndMoments.x()/k_.x();
    bodyState_[1] = forcesAndMoments.y()/k_.y();
    bodyState_[2] = forcesAndMoments.z()/ktheta_;

    Info << "Bodystate: x = " << bodyState_[0] << ", y = " << bodyState_[1] << ", theta = " << bodyState_[2] << endl;
    return vector(bodyState_[0],bodyState_[1],bodyState_[2]);
}*/

const Field<vectorField>& ThreeDoFRigidBodyControlledFlap::getPosition() const{
    return currentPosition_;
}

void ThreeDoFRigidBodyControlledFlap::update()
{
    if(localTimeIndex() < mesh_.time().timeIndex()){
        bodyStateOldOld_ = bodyStateOld_; //Set old old bodyState
        bodyStateOld_ = bodyState_; //Set old bodyState
        bodyStateVelOld_ = bodyStateVel_; //Set old bodyState
        forAll(prevPosition_,ipatch){
            prevPosition_[ipatch] = currentPosition_[ipatch];
        }
		
		betaOldOld_ = betaOld_;	
        betaOld_ = beta_;
        flapUpdated_=false;
    }
    bodyMotionFunction::update();//This updates localTimeIndex
}

void ThreeDoFRigidBodyControlledFlap::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());    
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
