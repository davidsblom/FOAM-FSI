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
    RBFMotionSolverFiles

Description
    Radial basis function motion solver using rigid body motion

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    RBFMotionSolverRigid.C

\*---------------------------------------------------------------------------*/

#ifndef fsiInterface_H
#define fsiInterface_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "dimensionedScalar.H"
#include "dimensionedType.H"
#include "primitiveFields.H"
#include "volFields.H"
#include "aeroForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RBFMotionSolverFiles Declaration
\*---------------------------------------------------------------------------*/

class fsiInterface
{
    // Private data

		//- Coupling dictionary
		IOdictionary dict_;
		
		//- Fluid mesh
		const fvMesh& mesh_;
	
		//- fsi subiteractions counter
		label fsiIter_;
		
		//- Minimum Fsi Subiterations
		label minFsiIter_;
	
		//- Maximum Fsi Subiterations
		label maxFsiIter_;
	
		//- Minimum fsi residual
		scalar minResidual_;
	
		//- Upper bound relaxation factor
		scalar relaxFsiUpper_;

		//- Lower bound relaxation factor
		scalar relaxFsiLower_;

		//- Initial relaxation factor
		scalar relaxFsi0_;
		
		//- relaxation factor
		scalar relaxFsi_;
	
		//- Previous relaxation factor
		scalar relaxFsiPrev_;
	
		//- Fsi coupling scheme
		word couplingScheme_;
	
		//- Boolean if write to file
		bool writeToFile_;

		//- Boolean if relaxFsi_ should be reused from previous time step. Default is false
		bool reuseRelaxFsi_;
		
		//- forces on fluid patches per patch
		Field<vectorField> fluidPatchForces_;

		//- forces on fluid patches in single vectorfield
		vectorField fluidForces_;
		
		//- previous forces on fluid patches in single vectorfield
		vectorField fluidForcesPrev_;
		
		//- residual
		vectorField residual_;
		vectorField residualPrev_;
		scalar normRes_;
		
		//- Output file
		OFstream ofFsiProp_;

		//- Patch IDs of moving objects
		const labelList& movingPatchIDs_;

		//- forces object
		aeroForces forcesObj_;
	
	// Private Member Functions

        //- Disallow default bitwise copy construct
        fsiInterface(const fsiInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fsiInterface&);
        
        //- Calculate relaxation factor based on chosen scheme
        void setRelexationParamter();
        
        //- Set forces based on pressure, velocity, viscosity and density
        void setInterfaceForce(const dimensionedScalar& rhoFluid, const dimensionedScalar& nu,const volVectorField& U,const volScalarField& p);
        
        //- Fill fluidPatchForces_ with relaxed forces from patchForces_
		void syncPatchForcesAfterRelax();

        //- Relaxate the forces @(k+1) based on relaxFsi_
        void relax();
        
        //- Return normalised residual
        scalar normResidual();

public:

    // Constructors

        //- Construct from polyMesh
        fsiInterface
        (
            const fvMesh& mesh,
			const IOdictionary& dict,
			const labelList& movingPatchIDs
        );

    // Destructor

        virtual ~fsiInterface();

    // Member Functions
    
		//- Set initial force based on this parameters
		void initializeForce(const dimensionedScalar& rhoFluid, const dimensionedScalar& nu,const volVectorField& U,const volScalarField& p);
		
		//- Return reference to forces on patches
		const Field<vectorField>& getPatchForces() const;

		//- Return reference to forceField
		const vectorField& getForces() const;

		//- Get sum of forces
		vector getTotalForce() const;
		
		//- Reset relaxation parameters to initial (relaxFsi0_)->used for new time step
		void reset();

		//- Update fsiInterface parameters based on new force->needed at end of subiteration
		void update(const dimensionedScalar& rhoFluid, const dimensionedScalar& nu,const volVectorField& U,const volScalarField& p);
		
		//- Check if converged or not (can also be based on maximum iter)
		bool converged();
		
		//- Return index of current subiteration
		label iter() const;
		
		//- Write subiteration properties
		void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
