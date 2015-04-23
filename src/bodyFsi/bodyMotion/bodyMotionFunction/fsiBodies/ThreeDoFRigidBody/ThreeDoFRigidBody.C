
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
    ThreeDoFRigidBody

\*---------------------------------------------------------------------------*/
#include "ThreeDoFRigidBody.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(ThreeDoFRigidBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, ThreeDoFRigidBody, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
void ThreeDoFRigidBody::writeState(const scalar& time)
{
	//Body
    ofBody_ << time << "\t" << bodyState_[0] << "\t" << bodyState_[1] << "\t" << bodyState_[2]
            << "\t" << bodyState_[3] << "\t" << bodyState_[4] << "\t" << bodyState_[5]
            << "\t" << bodyForces_.x() << "\t" << bodyForces_.y() << "\t" << bodyMoments_.z() << endl;
}

void ThreeDoFRigidBody::checkInterpolationRK4force()
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

vector ThreeDoFRigidBody::interpolateRK4term(const scalar rkCoeff,const vector valueOld, const vector value) const
{
    vector termRK4 = vector::zero;
    switch(interpolationRK4force_)
    {
        case 1:
            termRK4 = value;
            break;
        case 2:
            termRK4 = 0.5*(value + valueOld);
            break;
        case 4:
            termRK4 = rkCoeff*value + (1-rkCoeff)*valueOld;
            break;
        default:
            termRK4 = rkCoeff*value + (1-rkCoeff)*valueOld;
            break;
        break;
    }

    return termRK4;
}

SquareMatrix<scalar> ThreeDoFRigidBody::getRK4Coeffs() const
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

scalarField ThreeDoFRigidBody::bdf3Coeffs() const
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

scalarField ThreeDoFRigidBody::bdf2Coeffs() const
{
    scalar deltaT = mesh_.time().deltaT().value();
    scalar deltaT0 = 0.0;
    if (mesh_.time().timeIndex() < 2)
    {
        Info << nl << "oneDoFRigidBody: Doing Euler time integration in first time-step" << nl <<endl;
        deltaT0 = GREAT;
    }
    else
    {
        deltaT0 = mesh_.time().deltaT0().value();
    }

    //Set coefficients based on deltaT and deltaT0
    scalarField coeffs(3,0.0);
    coeffs[0] = 1 + deltaT/(deltaT + deltaT0);
    coeffs[2] = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    coeffs[1] = -1*(coeffs[0] + coeffs[2]);

    return coeffs;
}


vector ThreeDoFRigidBody::interpolateRK4forces(const scalar rkCoeff) const
{
    return interpolateRK4term(rkCoeff,bodyForcesOld_,bodyForces_);
}

vector ThreeDoFRigidBody::interpolateRK4moments(const scalar rkCoeff) const
{
    return interpolateRK4term(rkCoeff,bodyMomentsOld_,bodyMoments_);
}

void ThreeDoFRigidBody::calculateForcesAndMoments(const Field<vectorField>& patchForces)
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

    bodyForces_ = totalForces;
    bodyMoments_ = totalMoments;
}

vector ThreeDoFRigidBody::calculatePosition(const Field<vectorField>& patchForces)
{
    const scalar& dt = mesh_.time().deltaT().value();

    calculateForcesAndMoments(patchForces);

    //Start by saving previous iteration
    bodyStatePrevIter_ = bodyState_;

    //Set values of x,y,theta,dxdt,dydt,dthetadt to old values initially (every subiteration of fsi also)
    scalar xbody = bodyStateOld_[0];
    scalar ybody = bodyStateOld_[1];
    scalar theta = bodyStateOld_[2];
    scalar dxbody = bodyStateOld_[3];
    scalar dybody = bodyStateOld_[4];
    scalar dtheta = bodyStateOld_[5];

    if(integrationScheme_=="BDF1")
    {
        bool converged = false;
        label bdf1iters = 0;
        scalar res = 1;
        do{
            bdf1iters += 1;
            scalarField bodyStatePrev = bodyState_;
            //6 equations (for the 6 state variables)
            simpleMatrix<scalar> C(6,0,0);
            //Equation 1: d(x)/dt - dxdt = 0;
            C[0][0] = 1.0/dt; C[0][3] = -1;
            //Equation 2: d(y)/dt - dydt = 0;
            C[1][1] = 1.0/dt; C[1][4] = -1;
            //Equation 3: d(theta)/dt - dthetadt = 0;
            C[2][2] = 1.0/dt; C[2][5] = -1;
            //Equation 4: m*d(dxdt)/dt + c*dxdt + k*x = Fx + A2*dthetadt^2 + A1*d(dtheta)/dt
            scalar A1=Lcg2rc_*Foam::sin(theta+thetaGeom_);
            scalar A2=Lcg2rc_*Foam::cos(theta+thetaGeom_);
            C[3][0] = k_.x()/m_;
            C[3][3] = c_.x()/m_ + 1.0/dt;
            C[3][5] = -A1*1.0/dt - A2*dtheta;
            //Equation 5: m*d(dydt)/dt + c*dydt + k*y = Fy + A1*dthetadt^2 - A2*d(dtheta)/dt
            C[4][1] = k_.y()/m_;
            C[4][4] = c_.y()/m_ + 1.0/dt;
            C[4][5] = -A1*dtheta + A2*1.0/dt;
            //Equation 6:
            scalar Irc = Icg_+m_*pow(Lcg2rc_,2);
            C[5][2] = ktheta_/Irc;
            C[5][5] = ctheta_/Irc + 1.0/dt;
            C[5][3] = -m_/Irc*A1*1.0/dt;
            C[5][4] = m_/Irc*A2*1.0/dt;

            //RHS
            C.source()[0] = bodyStateOld_[0]*1.0/dt;
            C.source()[1] = bodyStateOld_[1]*1.0/dt;
            C.source()[2] = bodyStateOld_[2]*1.0/dt;
            C.source()[3] = bodyForces_.x()/m_ + bodyStateOld_[3]*1.0/dt - A1*bodyStateOld_[5]*1.0/dt;
            C.source()[4] = bodyForces_.y()/m_ + bodyStateOld_[4]*1.0/dt + A2*bodyStateOld_[5]*1.0/dt;
            C.source()[5] = bodyMoments_.z()/Irc + bodyStateOld_[5]*1.0/dt - m_/Irc*A1*bodyStateOld_[3]*1.0/dt + m_/Irc*A2*bodyStateOld_[4]*1.0/dt;

            //Solve system and update variables
            bodyState_ = C.LUsolve();
            xbody = bodyState_[0];
            ybody = bodyState_[1];
            theta = bodyState_[2];
            dxbody = bodyState_[3];
            dybody = bodyState_[4];
            dtheta = bodyState_[5];

            //Check for convergence
            res= ::sqrt(sum(sqr(bodyState_-bodyStatePrev)))/(::sqrt(sum(sqr(bodyState_)))+SMALL);
            converged = res<1e-10;
            //Info << "body BDF1: iter " <<bdf1iters <<", residual " << res << endl;
        }while(!converged);
        //Info << "body BDF1: converged in " << bdf1iters << " iterations. Final residual = " << res << endl;
    }
    else if(integrationScheme_=="BDF2")
    {
        //bodyState =[x,y,theta,dxdt,dydt,dthetadt]
        scalarField coeffs = bdf2Coeffs();

        bool converged = false;
        label bdf2iters = 0;
        scalar res = 1;
        do{
            bdf2iters += 1;
            scalarField bodyStatePrev = bodyState_;
            //6 equations (for the 6 state variables)
            simpleMatrix<scalar> C(6,0,0);
            //Equation 1: d(x)/dt - dxdt = 0;
            C[0][0] = coeffs[0]/dt; C[0][3] = -1;
            //Equation 2: d(y)/dt - dydt = 0;
            C[1][1] = coeffs[0]/dt; C[1][4] = -1;
            //Equation 3: d(theta)/dt - dthetadt = 0;
            C[2][2] = coeffs[0]/dt; C[2][5] = -1;
            //Equation 4: m*d(dxdt)/dt + c*dxdt + k*x = Fx + A2*dthetadt^2 + A1*d(dtheta)/dt
            scalar A1=Lcg2rc_*Foam::sin(theta+thetaGeom_);
            scalar A2=Lcg2rc_*Foam::cos(theta+thetaGeom_);
            C[3][0] = k_.x()/m_;
            C[3][3] = c_.x()/m_ + coeffs[0]/dt;
            C[3][5] = -A1*coeffs[0]/dt - A2*dtheta;
            //Equation 5: m*d(dydt)/dt + c*dydt + k*y = Fy + A1*dthetadt^2 - A2*d(dtheta)/dt
            C[4][1] = k_.y()/m_;
            C[4][4] = c_.y()/m_ + coeffs[0]/dt;
            C[4][5] = -A1*dtheta + A2*coeffs[0]/dt;
            //Equation 6:
            scalar Irc = Icg_+m_*pow(Lcg2rc_,2);
            C[5][2] = ktheta_/Irc;
            C[5][5] = ctheta_/Irc + coeffs[0]/dt;
            C[5][3] = -m_/Irc*A1*coeffs[0]/dt;
            C[5][4] = m_/Irc*A2*coeffs[0]/dt;

            //RHS
            C.source()[0] = -bodyStateOld_[0]*coeffs[1]/dt - bodyStateOldOld_[0]*coeffs[2]/dt;
            C.source()[1] = -bodyStateOld_[1]*coeffs[1]/dt - bodyStateOldOld_[1]*coeffs[2]/dt;
            C.source()[2] = -bodyStateOld_[2]*coeffs[1]/dt - bodyStateOldOld_[2]*coeffs[2]/dt;
            C.source()[3] = bodyForces_.x()/m_ - bodyStateOld_[3]*coeffs[1]/dt - bodyStateOldOld_[3]*coeffs[2]/dt + A1*(bodyStateOld_[5]*coeffs[1]/dt + bodyStateOldOld_[5]*coeffs[2]/dt);
            C.source()[4] = bodyForces_.y()/m_ - bodyStateOld_[4]*coeffs[1]/dt - bodyStateOldOld_[4]*coeffs[2]/dt - A2*(bodyStateOld_[5]*coeffs[1]/dt + bodyStateOldOld_[5]*coeffs[2]/dt);
            C.source()[5] = bodyMoments_.z()/Irc - bodyStateOld_[5]*coeffs[1]/dt - bodyStateOldOld_[5]*coeffs[2]/dt + m_/Irc*A1*(bodyStateOld_[3]*coeffs[1]/dt + bodyStateOldOld_[3]*coeffs[2]/dt)
                                - m_/Irc*A2*(bodyStateOld_[4]*coeffs[1]/dt + bodyStateOldOld_[4]*coeffs[2]/dt);

            //Solve system and update variables
            bodyState_ = C.LUsolve();
            xbody = bodyState_[0];
            ybody = bodyState_[1];
            theta = bodyState_[2];
            dxbody = bodyState_[3];
            dybody = bodyState_[4];
            dtheta = bodyState_[5];

            //Check for convergence
            res= ::sqrt(sum(sqr(bodyState_-bodyStatePrev)))/(::sqrt(sum(sqr(bodyState_)))+SMALL);
            converged = res<1e-10;
            //Info << "body BDF2: iter " <<bdf2iters <<", residual " << res << endl;
        }while(!converged);
        //Info << "body BDF2: converged in " << bdf2iters << " iterations. Final residual = " << res << endl;
    }
    else if(integrationScheme_=="BDF3")
    {
        //bodyState =[x,y,theta,dxdt,dydt,dthetadt]
        scalarField coeffs = bdf3Coeffs();

        bool converged = false;
        label bdf3iters = 0;
        scalar res = 1;
        do{
            //Set variables for subiterations
            bdf3iters += 1;
            scalarField bodyStatePrev = bodyState_;

            //6 equations (for the 6 state variables)
            simpleMatrix<scalar> C(6,0,0);
            //Equation 1: d(x)/dt - dxdt = 0;
            C[0][0] = coeffs[0]/dt; C[0][3] = -1;
            //Equation 2: d(y)/dt - dydt = 0;
            C[1][1] = coeffs[0]/dt; C[1][4] = -1;
            //Equation 3: d(theta)/dt - dthetadt = 0;
            C[2][2] = coeffs[0]/dt; C[2][5] = -1;
            //Equation 4: m*d(dxdt)/dt + c*dxdt + k*x = Fx + A2*dthetadt^2 + A1*d(dtheta)/dt
            scalar A1=Lcg2rc_*Foam::sin(theta+thetaGeom_);
            scalar A2=Lcg2rc_*Foam::cos(theta+thetaGeom_);
            C[3][0] = k_.x()/m_;
            C[3][3] = c_.x()/m_ + coeffs[0]/dt;
            C[3][5] = -A1*coeffs[0]/dt - A2*dtheta;
            //Equation 5: m*d(dydt)/dt + c*dydt + k*y = Fy + A1*dthetadt^2 - A2*d(dtheta)/dt
            C[4][1] = k_.y()/m_;
            C[4][4] = c_.y()/m_ + coeffs[0]/dt;
            C[4][5] = -A1*dtheta + A2*coeffs[0]/dt;
            //Equation 6:
            scalar Irc = Icg_+m_*pow(Lcg2rc_,2);
            C[5][2] = ktheta_/Irc;
            C[5][5] = ctheta_/Irc + coeffs[0]/dt;
            C[5][3] = -m_/Irc*A1*coeffs[0]/dt;
            C[5][4] = m_/Irc*A2*coeffs[0]/dt;

            //RHS
            C.source()[0] = -bodyStateOld_[0]*coeffs[1]/dt - bodyStateOldOld_[0]*coeffs[2]/dt - bodyStateOldOldOld_[0]*coeffs[3]/dt;
            C.source()[1] = -bodyStateOld_[1]*coeffs[1]/dt - bodyStateOldOld_[1]*coeffs[2]/dt - bodyStateOldOldOld_[1]*coeffs[3]/dt;;
            C.source()[2] = -bodyStateOld_[2]*coeffs[1]/dt - bodyStateOldOld_[2]*coeffs[2]/dt - bodyStateOldOldOld_[2]*coeffs[3]/dt;;
            C.source()[3] = bodyForces_.x()/m_ - bodyStateOld_[3]*coeffs[1]/dt - bodyStateOldOld_[3]*coeffs[2]/dt - bodyStateOldOldOld_[3]*coeffs[3]/dt
                            + A1*(bodyStateOld_[5]*coeffs[1]/dt + bodyStateOldOld_[5]*coeffs[2]/dt + bodyStateOldOldOld_[5]*coeffs[3]/dt);
            C.source()[4] = bodyForces_.y()/m_ - bodyStateOld_[4]*coeffs[1]/dt - bodyStateOldOld_[4]*coeffs[2]/dt - bodyStateOldOldOld_[4]*coeffs[3]/dt
                            - A2*(bodyStateOld_[5]*coeffs[1]/dt + bodyStateOldOld_[5]*coeffs[2]/dt + bodyStateOldOldOld_[5]*coeffs[3]/dt);
            C.source()[5] = bodyMoments_.z()/Irc - bodyStateOld_[5]*coeffs[1]/dt - bodyStateOldOld_[5]*coeffs[2]/dt - bodyStateOldOldOld_[5]*coeffs[3]/dt
                            + m_/Irc*A1*(bodyStateOld_[3]*coeffs[1]/dt + bodyStateOldOld_[3]*coeffs[2]/dt + bodyStateOldOldOld_[3]*coeffs[3]/dt)
                            - m_/Irc*A2*(bodyStateOld_[4]*coeffs[1]/dt + bodyStateOldOld_[4]*coeffs[2]/dt + bodyStateOldOldOld_[4]*coeffs[3]/dt);

            //Solve system and update variables
            bodyState_ = C.LUsolve();
            xbody = bodyState_[0];
            ybody = bodyState_[1];
            theta = bodyState_[2];
            dxbody = bodyState_[3];
            dybody = bodyState_[4];
            dtheta = bodyState_[5];

            //Check for convergence
            res= ::sqrt(sum(sqr(bodyState_-bodyStatePrev)))/(::sqrt(sum(sqr(bodyState_)))+SMALL);
            converged = res<1e-10;
            //Info << "body BDF2: iter " <<bdf2iters <<", residual " << res << endl;
        }while(!converged);
        //Info << "body BDF3: converged in " << bdf3iters << " iterations. Final residual = " << res << endl;
    }
    else if(integrationScheme_=="RK4")
    {
        SquareMatrix<scalar> rkMatrix = getRK4Coeffs();

        label nStages=rkMatrix.n()-1;
        Field<scalarField> bodyStateVelStages(rkMatrix.n()-1,scalarField(bodyState_.size(),0.0));
        for(int i=0;i<nStages;i++)
        {
            //Determine current state based on bodyStateVel
            scalarField bodyStageState(bodyStateOld_);
            for(int k=0;k<i;k++){
                bodyStageState += dt*rkMatrix[i][k+1]*bodyStateVelStages[k];
            }
            xbody=bodyStageState[0];ybody=bodyStageState[1];theta=bodyStageState[2];
            dxbody=bodyStageState[3];dybody=bodyStageState[4];dtheta=bodyStageState[5];

            //Interpolate forces
            vector bodyForces = interpolateRK4forces(rkMatrix[i][0]);
            vector bodyMoments = interpolateRK4moments(rkMatrix[i][0]);

            simpleMatrix<scalar> C(3,0,0);
            C[0][0] = m_;
            C[0][2] = -m_*Lcg2rc_*sin(theta+thetaGeom_);
            C[1][1] = m_;
            C[1][2] = m_*Lcg2rc_*cos(theta+thetaGeom_);
            C[2][0] = -m_*Lcg2rc_*sin(theta+thetaGeom_);
            C[2][1] = m_*Lcg2rc_*cos(theta+thetaGeom_);
            C[2][2] = Icg_ + m_*pow(Lcg2rc_,2);

            C.source()[0] = bodyForces.x() - c_.x()*dxbody - k_.x()*xbody + m_*Lcg2rc_*pow(dtheta,2)*cos(theta+thetaGeom_);
            C.source()[1] = bodyForces.y() - c_.y()*dybody - k_.y()*ybody + m_*Lcg2rc_*pow(dtheta,2)*sin(theta+thetaGeom_);
            C.source()[2] = bodyMoments.z() - ctheta_*dtheta - ktheta_*theta;//Add some moment here
            List<scalar> bodyStateAcc = C.LUsolve();

            //Calculate bodyStateVel for this stage.
            bodyStateVelStages[i][0]=dxbody;
            bodyStateVelStages[i][1]=dybody;
            bodyStateVelStages[i][2]=dtheta;
            bodyStateVelStages[i][3]=bodyStateAcc[0];
            bodyStateVelStages[i][4]=bodyStateAcc[1];
            bodyStateVelStages[i][5]=bodyStateAcc[2];
        }
        //calculate final state at t(n+1)
        bodyState_ = bodyStateOld_;
        for(int i=0;i<nStages;i++){
            bodyState_ += dt*bodyStateVelStages[i]*rkMatrix[4][i+1];
        }
    }
    else
    {
        FatalErrorIn("void ThreeDoFRigidBody::calculatePosition()")
            << "integrationScheme " << integrationScheme_ << " not available."
            << "Available schemes are: (RK4, BDF1, BDF2, BDF3)"
            << abort(FatalError);
    }

    vector bodyStatePos(bodyState_[0],bodyState_[1],bodyState_[2]);
    Info << "bodyState = " << bodyStatePos << endl;

    return bodyStatePos;
}

void ThreeDoFRigidBody::readInitialStateIfPresent()
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

    bodyStateOldOldOld_ = bodyStateOldOld_;
    bodyStateOldOld_ = bodyStateOld_;
    bodyStateOld_ = bodyState_;

    bodyState0_ = bodyState_;

    Info << "Initial state is set to: " << bodyState_ << endl;
}

void ThreeDoFRigidBody::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void ThreeDoFRigidBody::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}

// * * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

tmp<Field<vectorField> > ThreeDoFRigidBody::getChildMotion()
{
    FatalErrorIn("void ThreeDoFRigidBody::getChildMotion()")
            << "this bodyMotionFunction cannot be used in a non-fsi solver, since getChildMotion(const Field<vectorField>& patchForces) needs to be called"
            << abort(FatalError);

    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    //Field<vectorField>& motion = tMotion();

    return tMotion;
}

tmp<Field<vectorField> > ThreeDoFRigidBody::getChildMotion(const Field<vectorField>& patchForces)
{
    Info << "Solve Solid " << name_ << endl;
    //Prepare return variable
    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    //Calculate the bodyState based on the forces
    calculatePosition(patchForces);

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

    return tMotion;
}

tmp<Field<vectorField> > ThreeDoFRigidBody::getChildSteadyStateMotion(const Field<vectorField>& patchForces)
{
     Info << "Solve Solid " << name_ << endl;
    //Prepare return variable
    tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    calculateForcesAndMoments(patchForces);

    //Start by saving previous iteration
    bodyStatePrevIter_ = bodyState_;

    //Calculate steady state deformation
    bodyState_[0] = bodyForces_.x()/k_.x();
    bodyState_[1] = bodyForces_.y()/k_.y();
    bodyState_[2] = bodyMoments_.z()/ktheta_;

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

    return tMotion;
}

void ThreeDoFRigidBody::setInitialPoints(){
    forAll(patchIDs_,ipatch){
        initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
    }
}

const labelList ThreeDoFRigidBody::getChildMotionPatchIDs() const{
    return patchIDs_;
}

const wordList ThreeDoFRigidBody::getChildMotionPatchNames() const{
    return patchNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ThreeDoFRigidBody::ThreeDoFRigidBody(const dictionary& dict,const fvMesh& mesh,const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
integrationScheme_(dict_.lookup("integrationScheme")),
interpolationRK4force_(dict_.lookupOrDefault("interpolation",4)),
k_(dict_.lookup("k")),
c_(dict_.lookup("c")),
ktheta_(readScalar(dict_.lookup("ktheta"))),
ctheta_(readScalar(dict_.lookup("ctheta"))),
m_(0.0),
Icg_(readScalar(dict_.lookup("Icg"))),
Lcg2rc_(readScalar(dict_.lookup("Lcg2rc"))),
rc_(dict_.lookup("rc")),
thetaGeom_(readScalar(dict_.lookup("thetaGeom"))),
bodyState_(6,0.0),
bodyStatePrevIter_(6,0.0),
bodyStateOld_(6,0.0),
bodyStateOldOld_(6,0.0),
bodyStateOldOldOld_(6,0.0),
bodyState0_(6,0.0),
theta0_(0),
bodyForces_(vector::zero),
bodyForcesOld_(vector::zero),
bodyMoments_(vector::zero),
bodyMomentsOld_(vector::zero),
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

    //Read initial state
    readInitialStateIfPresent();

    if(writeToFile())
    {
        ofBody_ << "time\tx\ty\ttheta\tdxdt\tdydt\tdthetadt\tFx\tFy\tMtheta" << endl;
    }
    else
    {
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ThreeDoFRigidBody::~ThreeDoFRigidBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& ThreeDoFRigidBody::getPosition() const{
    return currentPosition_;
}

void ThreeDoFRigidBody::update()
{
    if(localTimeIndex() < mesh_.time().timeIndex()){
        bodyForcesOld_ = bodyForces_;
        bodyMomentsOld_ = bodyMoments_;
        bodyStateOldOldOld_ = bodyStateOldOld_; //Set old old old bodyState
        bodyStateOldOld_ = bodyStateOld_; //Set old old bodyState
        bodyStateOld_ = bodyState_; //Set old bodyState
        forAll(prevPosition_,ipatch){
            prevPosition_[ipatch] = currentPosition_[ipatch];
        }
    }
    bodyMotionFunction::update();//This updates localTimeIndex
}

void ThreeDoFRigidBody::write()
{
    if(writeToFile())
    {
        writeState(mesh_.time().value());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
