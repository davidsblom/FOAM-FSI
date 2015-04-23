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
    oneDoFRigidBody

\*---------------------------------------------------------------------------*/
#include "oneDoFRigidBody.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(oneDoFRigidBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, oneDoFRigidBody, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
void oneDoFRigidBody::writeState(const scalar& time)
{
    //Body
    ofBody_ << time << "\t" << bodyState_[0] << "\t" << bodyState_[1] << "\t" << bodyForce_ << endl;
}

void oneDoFRigidBody::checkInterpolationRK4force()
{
    if(integrationScheme_=="RK4"){
        switch(interpolationRK4force_)
        {
            case 1:
                Info << "oneDoFRigidBody: Doing force interpolation at t(n+1)" << endl;
                break;
            case 2:
                Info << "oneDoFRigidBody: Doing force interpolation at 0.5*[t(n+1)+t(n)]" << endl;
                break;
            case 4:
                Info << "oneDoFRigidBody: Doing force interpolation at [k*t(n+1) + (1-k)*t(n)]" << endl;
                break;
            default:
                Info << "oneDoFRigidBody: Unknown force interpolation for RK4. Setting to best." << endl;
                interpolationRK4force_ = 4;
                break;
        }
    }
}

scalar oneDoFRigidBody::interpolateRK4force(const scalar rkCoeff) const
{
    scalar forceRK4 = 0.0;
    switch(interpolationRK4force_)
    {
        case 1:
            forceRK4 = bodyForce_;
            break;
        case 2:
            forceRK4 = 0.5*(bodyForce_ + bodyForceOld_);
            break;
        case 4:
            forceRK4 = rkCoeff*bodyForce_ + (1-rkCoeff)*bodyForceOld_;
            break;
        default:
            forceRK4 = rkCoeff*bodyForce_ + (1-rkCoeff)*bodyForceOld_;
            break;
        break;
    }

    return forceRK4;
}

SquareMatrix<scalar> oneDoFRigidBody::getRK4Coeffs() const
{
    // === matrix 5x5 === //
    // 0     0     0     0     0
    // 1/3   1/3   0     0     0
    // 2/3   -1/3  1     0     0
    // 1     1     -1    1     0
    // 0     1/8   3/8   3/8   1/8
    // === end === //
    SquareMatrix<scalar> rkMatrix(5,0.0);
    rkMatrix[1][0]=1.0/3.0;rkMatrix[2][0]=2.0/3.0;rkMatrix[3][0]=1.0;
    rkMatrix[4][1]=1.0/8.0;rkMatrix[4][2]=3.0/8.0;rkMatrix[4][3]=3.0/8.0;rkMatrix[4][4]=1.0/8.0;
    rkMatrix[1][1]=1.0/3.0;
    rkMatrix[2][1]=-1.0/3.0;rkMatrix[2][2]=1.0;
    rkMatrix[3][1]=1.0;rkMatrix[3][2]=-1.0;rkMatrix[3][3]=1.0;

    /*// === matrix 5x5 === //
    // 0     0     0     0     0
    // 1/2   1/2   0     0     0
    // 1/2   0     1/2   0     0
    // 1     0     0     1     0
    // 0     1/6   1/3   1/3   1/6
    // === end === //
    SquareMatrix<scalar> rkMatrix(5,0.0);
    rkMatrix[1][0]=1.0/2.0;rkMatrix[2][0]=1.0/2.0;rkMatrix[3][0]=1.0;
    rkMatrix[4][1]=1.0/6.0;rkMatrix[4][2]=1.0/3.0;rkMatrix[4][3]=1.0/3.0;rkMatrix[4][4]=1.0/6.0;
    rkMatrix[1][1]=1.0/2.0;
    rkMatrix[2][2]=1.0/2.0;
    rkMatrix[3][3]=1.0;*/

    return rkMatrix;
}

scalarField oneDoFRigidBody::bdf3Coeffs() const
{
    scalarField coeffs(4,0.0);

    //if first time step
    if (mesh_.time().timeIndex() < 2)
    {
        coeffs[0]=1.0;
        coeffs[1]=-1.0;
    }
    else if(mesh_.time().timeIndex() < 3)
    {
        coeffs[0]=3.0/2.0;
        coeffs[1]=-2.0;
        coeffs[2]=1.0/2.0;
    }
    else
    {
        coeffs[0]=11.0/6.0;
        coeffs[1]=-3.0;
        coeffs[2]=3.0/2.0;
        coeffs[3]=-1.0/3.0;
    }

    return coeffs;
}

scalar oneDoFRigidBody::deltaT_() const
{
    return mesh_.time().deltaT().value();
}


scalar oneDoFRigidBody::deltaT0_() const
{
    if (mesh_.time().timeIndex() < 2)
    {
        Info << nl << "oneDoFRigidBody: Doing Euler time integration in first time-step" << nl <<endl;
        return GREAT;
    }
    else
    {
        return mesh_.time().deltaT0().value();
    }
}

scalar oneDoFRigidBody::calculateForce(const Field<vectorField>& patchForces) const
{
    scalar force(0);

    vector totalForces(vector::zero);
    forAll(patchIDs_,ipatch)
    {
        totalForces += gSum(patchForces[patchIDs_[ipatch]]);
    }

    force = totalForces.y();

    Info << "force = " << force << endl;

    return force;
}

scalar oneDoFRigidBody::calculatePosition(const Field<vectorField>& patchForces)
{
    Info << "Solve Solid" << endl;
    const scalar& dt = mesh_.time().deltaT().value();

    bodyForce_ = calculateForce(patchForces);

    //Start by saving previous iteration
    bodyStatePrevIter_ = bodyState_;

    //Set values of x,y,theta,dxdt,dydt,dthetadt to old values initially (every subiteration of fsi also)
    scalar ybody = bodyStateOld_[0];
    scalar dybody = bodyStateOld_[1];

    if(integrationScheme_=="RK4")
    {
        SquareMatrix<scalar> rkMatrix = getRK4Coeffs();

        label nStages=rkMatrix.n()-1;
        Field<scalarField> bodyStateVelStages(rkMatrix.n()-1,scalarField(2,0.0));
        for(int i=0;i<nStages;i++)
        {
            //Determine current state based on bodyStateVel
            scalarField bodyStageState(bodyStateOld_);
            for(int k=0;k<i;k++){
                bodyStageState += dt*rkMatrix[i][k+1]*bodyStateVelStages[k];
            }
            ybody=bodyStageState[0];dybody=bodyStageState[1];

            //Interpolate forces
            //scalar bodyForce = bodyForce_;
            //scalar bodyForce = 0.5*(bodyForce_ + bodyForceOld_);
            //scalar bodyForce = rkMatrix[i][0]*bodyForce_ + (1-rkMatrix[i][0])*bodyForceOld_;
            scalar bodyForce = interpolateRK4force(rkMatrix[i][0]);

            //Calculate bodyStateVel for this stage.
            bodyStateVelStages[i][0]=dybody;
            bodyStateVelStages[i][1]=(bodyForce/m_ - c_/m_*dybody - k_/m_*ybody);
            Info << endl;
        }
        bodyState_ = bodyStateOld_;
        for(int i=0;i<nStages;i++){
            bodyState_ += dt*bodyStateVelStages[i]*rkMatrix[4][i+1];
        }
    }
    else if(integrationScheme_=="BDF1")
    {
        //Set old body state variables
        scalar ybody0 = bodyStateOld_[0];
        scalar dybody0 = bodyStateOld_[1];

        //Matrix with equations
        simpleMatrix<scalar> A(2,0,0);
        A[0][0] = c_/m_;
        A[0][1] = k_/m_;
        A[1][0] = -1;
        A[1][1] = 0;

        //Matrix with diagonal BDF1 terms
        A[0][0] += 1/dt;
        A[1][1] += 1/dt;

        //rhs
        A.source()[0] = bodyForce_/m_ - (-1/dt*dybody0);
        A.source()[1] = - (-1/dt*ybody0);

        //solve, where solution is dybody,ybody
        List<scalar> bodyState = A.LUsolve();
        bodyState_[1] = bodyState[0];
        bodyState_[0] = bodyState[1];
    }
    else if(integrationScheme_=="BDF2")
    {
        scalar deltaT = deltaT_();
        scalar deltaT0 = deltaT0_();
        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        scalar coefft0  = -1*(coefft + coefft00);

        //Set old body state variables
        scalar ybody0 = bodyStateOld_[0];
        scalar dybody0 = bodyStateOld_[1];
        scalar ybody00 = bodyStateOldOld_[0];
        scalar dybody00 = bodyStateOldOld_[1];

        //Matrix with equations
        simpleMatrix<scalar> A(2,0,0);
        A[0][0] = c_/m_;
        A[0][1] = k_/m_;
        A[1][0] = -1;
        A[1][1] = 0;

        //Matrix with diagonal BDF2 terms
        A[0][0] += coefft/dt;
        A[1][1] += coefft/dt;

        //rhs
        A.source()[0] = bodyForce_/m_ - (coefft0/dt*dybody0 + coefft00*dybody00/dt);
        A.source()[1] = - (coefft0/dt*ybody0 + coefft00*ybody00/dt);

        //solve, where solution is dybody,ybody
        List<scalar> bodyState = A.LUsolve();
        bodyState_[1] = bodyState[0];
        bodyState_[0] = bodyState[1];
    }
    else if(integrationScheme_=="BDF3")
    {
        scalarField coeffs = bdf3Coeffs();

        //Set old body state variables
        scalar ybody0 = bodyStateOld_[0];
        scalar dybody0 = bodyStateOld_[1];
        scalar ybody00 = bodyStateOldOld_[0];
        scalar dybody00 = bodyStateOldOld_[1];
        scalar ybody000 = bodyStateOldOldOld_[0];
        scalar dybody000 = bodyStateOldOldOld_[1];

        //Matrix with equations
        simpleMatrix<scalar> A(2,0,0);
        A[0][0] = c_/m_;
        A[0][1] = k_/m_;
        A[1][0] = -1;
        A[1][1] = 0;

        //Matrix with diagonal BDF2 terms
        A[0][0] += coeffs[0]/dt;
        A[1][1] += coeffs[0]/dt;

        //rhs
        A.source()[0] = bodyForce_/m_ - (coeffs[1]/dt*dybody0 + coeffs[2]*dybody00/dt + coeffs[3]*dybody000/dt);
        A.source()[1] = - (coeffs[1]/dt*ybody0 + coeffs[2]*ybody00/dt + coeffs[3]*ybody000/dt);

        //solve, where solution is dybody,ybody
        List<scalar> bodyState = A.LUsolve();
        bodyState_[1] = bodyState[0];
        bodyState_[0] = bodyState[1];
    }
    else
    {
        FatalErrorIn("void oneDoFRigidBody::calculatePosition()")
            << "integrationScheme " << integrationScheme_ << " not available."
            << "Available schemes are: (RK4,BDF1, BDF2, BDF3)"
            << abort(FatalError);
    }

    return bodyState_[0];
}

void oneDoFRigidBody::readInitialStateIfPresent()
{
    scalar state(dict_.lookupOrDefault("initialState",0.0));
    scalar stateVel(dict_.lookupOrDefault("initialStateVelocity",0.0));

    bodyState_[0] = state;
    bodyState_[1] = stateVel;

    bodyStateOld_ = bodyState_;
    bodyStateOldOld_ = bodyStateOld_;
    bodyStateOldOldOld_ = bodyStateOldOld_;

    bodyState0_ = bodyState_;

    Info << "Initial state is set to: " << bodyState_ << endl;
}

void oneDoFRigidBody::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void oneDoFRigidBody::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

// * * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

tmp<Field<vectorField> > oneDoFRigidBody::getChildMotion()
{
    FatalErrorIn("void oneDoFRigidBody::getChildMotion()")
            << "this bodyMotionFunction cannot be used in a non-fsi solver, since getChildMotion(const Field<vectorField>& patchForces) needs to be called"
            << abort(FatalError);

    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));

    return tMotion;
}

tmp<Field<vectorField> > oneDoFRigidBody::getChildMotion(const Field<vectorField>& patchForces)
{
    //Prepare return variable
    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    //Calculate the bodyState based on the forces
    calculatePosition(patchForces);

    //Translation is subtracting current state with previous state
    vector translationVector(0,bodyState_[0] - bodyStatePrevIter_[0],0);

    forAll(patchIDs_,ipatch){
        motion[ipatch] = vectorField(initialPoints_[ipatch].size(),vector::zero);
        motion[ipatch] = translationVector;
    }



    return tMotion;
}

void oneDoFRigidBody::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
    }
}

const labelList oneDoFRigidBody::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList oneDoFRigidBody::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneDoFRigidBody::oneDoFRigidBody(const dictionary& dict,const fvMesh& mesh,const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
integrationScheme_(dict_.lookup("integrationScheme")),
interpolationRK4force_(dict_.lookupOrDefault("interpolation",4)),
k_(readScalar(dict_.lookup("k"))),
c_(readScalar(dict_.lookup("c"))),
m_(0.0),
bodyState_(2,0.0),
bodyStatePrevIter_(2,0.0),
bodyStateOld_(2,0.0),
bodyStateOldOld_(2,0.0),
bodyStateOldOldOld_(2,0.0),
bodyState0_(2,0.0),
bodyForce_(0),
bodyForceOld_(0),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
    setPatchIDs();
    setInitialPoints();
    checkInterpolationRK4force();

    //Determine thickness of mesh and set m accordingly
    Vector<label> e = (-mesh.geometricD()+Vector<label>(1,1,1))/2;
    scalar meshDepth(vector(e.x(),e.y(),e.z()) & vector(mesh.bounds().span()));
    scalar m(readScalar(dict_.lookup("m")));
    m_ = m*meshDepth;

    Info << "BodyProperties: m = " << m_ << ", k = " << k_ <<", c = " << c_ << endl;

    //Read initial state
    readInitialStateIfPresent();

    if(writeToFile())
    {
        ofBody_ << "time\ty\tdydt\tFy" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oneDoFRigidBody::~oneDoFRigidBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*vector oneDoFRigidBody::calculateSteadyStatePosition(const vectorField& forces)
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

const Field<vectorField>& oneDoFRigidBody::getPosition() const{
    return currentPosition_;
}

void oneDoFRigidBody::update()
{
    if(localTimeIndex() < mesh_.time().timeIndex()){
        bodyStateOldOldOld_ = bodyStateOldOld_;//set t(n-2) bodyState
        bodyStateOldOld_ = bodyStateOld_;//set t(n-1) bodyState
        bodyStateOld_ = bodyState_; //Set t(n) bodyState
        bodyForceOld_ = bodyForce_;
        forAll(prevPosition_,ipatch){
            prevPosition_[ipatch] = currentPosition_[ipatch];
        }
    }
    bodyMotionFunction::update();//This updates localTimeIndex
}

void oneDoFRigidBody::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
