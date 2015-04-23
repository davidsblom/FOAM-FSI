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
    RBFMotionSolverExt

Description
    Radial basis function motion solver using rigid body motion

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    RBFMotionSolverExt.C

\*---------------------------------------------------------------------------*/

#ifndef RBFMotionSolverExt_H
#define RBFMotionSolverExt_H

#include "motionSolver.H"
#include "polyMesh.H"
#include "RBFInterpolationReduced.H"
#include "twoDPointCorrectorRBF.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RBFMotionSolverExt Declaration
\*---------------------------------------------------------------------------*/

class RBFMotionSolverExt
:
    public motionSolver
{
    // Private data

        //- Static patch names
        wordList staticPatches_;

        //- Stathc patch IDs
        labelList staticPatchIDs_;

        //- Static patch names
        wordList movingPatches_;

        //- Moving patch IDs
        labelList movingPatchIDs_;

        //- Coarsening ratio for static boundary
        label coarseningRatioStat_;
        
        //- Coarsening ratio for static boundary
        label coarseningRatioMov_;

        //- Control point IDs
        labelList controlIDs_;

        //- Control points on the boundary
        mutable vectorField controlPoints_;

        //- Internal point IDs
        labelList internalIDs_;

        //- Internal points
        vectorField internalPoints_;

        //- Static point IDs
        labelList staticIDs_;

        //- Static points
        vectorField staticPoints_;

        //- Moving point IDs
        labelList movingIDs_;

        //- Moving points original location
        vectorField movingPoints_;

        //- Moving boolean for all patches (empty labelField if not moving at all)
        Field<labelField> patchMovingGamma_;
        
        //- Motion of moving points
        vectorField motion_;
        
        //- Index of which of the moving points are used as control points
        labelList movingControlPointsIndex_;

        //- point IDs of the moving control points for each PROC
        labelList movingControlIDs_;

		//- Index of where new patch starts in movingPoints_
		labelList movingPointsPatchStart_;

		//- Number of static control points (nMovingControlPoints = controlIDs_.size() - nStaticControlPoints_)
		label nStaticControlPoints_;

        //- List of moving control points (list is size of no of procs)
        labelList nGlobalMovingControlPoints_;
        labelList nGlobalStaticControlPoints_;
		
        //RBFInterpolation interpolationRed_;
        RBFInterpolationReduced interpolationRed_;

        //- New 2D correct thingie
        twoDPointCorrectorRBF twoDPointCorrectorRBF_;

// Private Member Functions

        //- Disallow default bitwise copy construct
        RBFMotionSolverExt(const RBFMotionSolverExt&);

        //- Disallow default bitwise assignment
        void operator=(const RBFMotionSolverExt&);


        //- Make control point IDs.  Constructor helper
        void makeControlIDs();

        //- Make control point IDs for 2D
        void makeControlIDs2D();

        //- Parallel processer controlPoint gatherer
        void gatherControlPoints();

        //- Returns the new point locations in the correct format to the solver
		void moveMesh(pointField& newPoints) const;

        void setAndOrderPatchIDs();

        //- Returns list of point IDs of the moving control points
        void setMovingControlPointsIDs();

        //- Returns list of indicies of moving control points in controlPoints_
        const labelList movingControlPointsIndices() const;


        void gatherControlPointsExt();

        //- Prepare for parallel run
        void prepareParallelRun();

        void checkInterpolationDimension();

public:

    //- Runtime type information
    TypeName("RBFMotionSolverExt");

    // Constructors

        //- Construct from polyMesh
        RBFMotionSolverExt
        (
            const polyMesh& mesh,
            Istream& msData
        );

    // Destructor

        virtual ~RBFMotionSolverExt();

    // Member Functions

		//- Set motion_ according to m
		void setMotion(const vectorField& m);

        //- Set motion_ according to m, where m is equal to number of patches and with non-empty vectorFields for moving patches
        void setMotion(const Field<vectorField>& m);
		
		//- Get movingPoints
		const vectorField& movingPoints() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
