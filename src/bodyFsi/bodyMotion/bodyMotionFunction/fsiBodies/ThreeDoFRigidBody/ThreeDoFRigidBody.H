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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    ThreeDoFRigidBody

Description
	Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    ThreeDoFRigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef ThreeDoFRigidBody_H
#define ThreeDoFRigidBody_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"
#include "simpleMatrix.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ThreeDoFRigidBody Declaration
\*---------------------------------------------------------------------------*/

class ThreeDoFRigidBody
:
    public bodyMotionFunction
{
	// Private data
	const dictionary& dict_;
	const fvMesh& mesh_;
    const word name_;

    //Patches of this body
    wordList patchNames_;
    labelList patchIDs_;

    //Solve properties
    word integrationScheme_;
    label interpolationRK4force_;

	//Body properties
	vector k_;
    vector c_;
    scalar ktheta_;
    scalar ctheta_;
    scalar m_;
    scalar Icg_;
    scalar Lcg2rc_;
    vector rc_;
    scalar thetaGeom_;//In radians

	//Body movement parameters: bodyState = [x,y,theta,dxdt,dydt,dthetadt]
    scalarField bodyState_;
    scalarField bodyStatePrevIter_;
    scalarField bodyStateOld_;
    scalarField bodyStateOldOld_;
    scalarField bodyStateOldOldOld_;
    scalarField bodyState0_;
    scalar theta0_;

    //Body forces
    vector bodyForces_;
    vector bodyForcesOld_;
    vector bodyMoments_;
    vector bodyMomentsOld_;

    //Positions
    Field<vectorField> currentPosition_;
    Field<vectorField> prevPosition_;
    Field<vectorField> initialPoints_;

    //Writing properties
    OFstream ofBody_;

    // Private Member Functions

        //- Writing function
        void writeState(const scalar& time);

        //- Check if interpolation type is set correct
        void checkInterpolationRK4force();

        //- Actual interpolation of forces or moments
        vector interpolateRK4term(const scalar rkCoeff,const vector valueOld, const vector value) const;

        //- Return the Runge-Kutta 4 coefficient matrix
        SquareMatrix<scalar> getRK4Coeffs() const;

        //- Return list of BDF3 coefficients (different for first two time steps)
        scalarField bdf3Coeffs() const;

        //- Return list of BDF2 coefficients (different for first two time steps)
        scalarField bdf2Coeffs() const;

        //- Return the interpolated force
        vector interpolateRK4forces(const scalar rkCoeff) const;

        //- Return the interpolated moments
        vector interpolateRK4moments(const scalar rkCoeff) const;

        //- Calcaulte moment around current rotation center (rc_ + (x,y,rc_.z))
        void calculateForcesAndMoments(const Field<vectorField>& patchForces);

        //- Actual position calculation
        vector calculatePosition(const Field<vectorField>& patchForces);

        //- Read initial state from dict if present (state = [x,y,theta,dxdt,dydt,dthetadt])
        void readInitialStateIfPresent();

        //- Set patchIDs_ based on the patchNames
        void setPatchIDs();

    	//- Disallow default bitwise copy construct
        ThreeDoFRigidBody(const ThreeDoFRigidBody&);

    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        virtual tmp<Field<vectorField> > getChildMotion(const Field<vectorField>& patchForces);

        virtual tmp<Field<vectorField> > getChildSteadyStateMotion(const Field<vectorField>& patchForces);

        void setInitialPoints();

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;

	public:

    //- Runtime type information
    TypeName("ThreeDoFRigidBody");

	//Public data

    // Constructors

		//- Construct from dict + mesh + name
        ThreeDoFRigidBody(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new ThreeDoFRigidBody(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~ThreeDoFRigidBody();

    // Member Functions

        //- Calculate steady state position based on forces
        //virtual vector calculateSteadyStatePosition(const vectorField& forces);

        virtual const Field<vectorField>& getPosition() const;

        virtual void update();

        //- Write out body properties if requested
        virtual void write();
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
