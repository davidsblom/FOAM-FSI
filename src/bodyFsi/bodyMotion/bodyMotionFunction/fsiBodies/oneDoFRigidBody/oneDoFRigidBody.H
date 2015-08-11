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
    oneDoFRigidBody

Description
    Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    oneDoFRigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef oneDoFRigidBody_H
#define oneDoFRigidBody_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bodyMotionFunction.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"
#include "simpleMatrix.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oneDoFRigidBody Declaration
\*---------------------------------------------------------------------------*/

class oneDoFRigidBody
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
    scalar k_;
    scalar c_;
    scalar m_;

    //Body movement parameters: bodyState = [x,y,theta,dxdt,dydt,dthetadt]
    scalarField bodyState_;
    scalarField bodyStatePrevIter_;
    scalarField bodyStateOld_;
    scalarField bodyStateOldOld_;
    scalarField bodyStateOldOldOld_;
    scalarField bodyState0_;

    //Body forces
    scalar bodyForce_;
    scalar bodyForceOld_;

    //Positions
    Field<vectorField> currentPosition_;
    Field<vectorField> prevPosition_;
    Field<vectorField> initialPoints_;

    //Writing properties
    OFstream ofBody_;

    // Private Member Functions

        //- Writing function
        void writeState(const scalar& time);

        //- Check if interpolation option exists
        void checkInterpolationRK4force();

        //- Interpolate function of RK force
        scalar interpolateRK4force(const scalar rkCoeff) const;

        //- Return RK4 Coeffs
        SquareMatrix<scalar> getRK4Coeffs() const;

        //- Return BDF3 coeffs
        scalarField bdf3Coeffs() const;

        //- Return deltaT
        scalar deltaT_() const;

        //- Return deltaT0 or GREAT if timeIndex < 2 (first time step)
        scalar deltaT0_() const;

        //- Calcaulte moment around current rotation center (rc_ + (x,y,rc_.z))
        scalar calculateForce(const Field<vectorField>& patchForces) const;

        scalar calculatePosition(const Field<vectorField>& patchForces);

        //- Read initial state from dict if present (state = [x,y,theta,dxdt,dydt,dthetadt])
        void readInitialStateIfPresent();

        //- Set patchIDs_ based on the patchNames
        void setPatchIDs();

        //- Disallow default bitwise copy construct
        oneDoFRigidBody(const oneDoFRigidBody&);

    protected:

        virtual tmp<Field<vectorField> > getChildMotion();

        virtual tmp<Field<vectorField> > getChildMotion(const Field<vectorField>& patchForces);

        void setInitialPoints();

        virtual const labelList getChildMotionPatchIDs() const;

        virtual const wordList getChildMotionPatchNames() const;

    public:

    //- Runtime type information
    TypeName("oneDoFRigidBody");

    //Public data

    // Constructors

        //- Construct from dict + mesh + name
        oneDoFRigidBody(const dictionary& dict,const fvMesh& mesh,const word name);

        virtual autoPtr<bodyMotionFunction> clone() const
        {
            return autoPtr<bodyMotionFunction>(new oneDoFRigidBody(this->dict_,this->mesh_,this->name_));
        };

    // Destructor

        virtual ~oneDoFRigidBody();

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
