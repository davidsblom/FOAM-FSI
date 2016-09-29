
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */


#ifndef foamFluidSolver_H
#define foamFluidSolver_H

#include <memory>

#include "BaseMultiLevelSolver.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "RBFMeshMotionSolver.H"
#include "twoDPointCorrectorRBF.H"
#include <time.h>
#include <unordered_map>
#include <vector>

using namespace fsi;
using std::shared_ptr;

class foamFluidSolver : public BaseMultiLevelSolver
{
    public:
        foamFluidSolver(
            const std::string & name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        ~foamFluidSolver();

        virtual void getAcousticsDensityLocal( matrix & data ) = 0;

        virtual void getAcousticsTemperatureLocal( matrix & data ) = 0;

        virtual void getAcousticsVelocityLocal( matrix & data ) = 0;

        virtual void getAcousticsPressureLocal( matrix & data ) = 0;

        virtual void getAcousticsTemperatureGradientLocal( matrix & data ) = 0;

        virtual void getAcousticsVelocityGradientLocal( matrix & data ) = 0;

        virtual void getAcousticsPressureGradientLocal( matrix & data ) = 0;

        virtual void getTractionLocal( matrix & traction ) = 0;

        virtual void getWritePositionsLocalAcoustics( matrix & writePositions ) = 0;

        virtual void finalizeTimeStep();

        virtual void getReadPositions( matrix & readPositions );

        virtual void getWritePositions( matrix & writePositions );

        virtual void moveMesh();

        virtual void solve() = 0;

        virtual void setAcousticsTemperatureGradient( const matrix & temperatureGradient ) = 0;

        virtual void setAcousticsPressureGradient( const matrix & pressureGradient ) = 0;

        virtual void setAcousticsVelocityGradient( const matrix & velocityGradient ) = 0;

        virtual void setAcousticsTemperature( const matrix & temperature ) = 0;

        virtual void setAcousticsPressure( const matrix & pressure ) = 0;

        virtual void setAcousticsVelocity( const matrix & velocity ) = 0;

        virtual void solve(
            const matrix & input,
            matrix & output
            );

        void setDisplacementLocal( const matrix & displacement );

        void getReadPositionsLocal( matrix & readPositions );

        void getWritePositionsLocal( matrix & writePositions );

        void run();

        std::string name;
        shared_ptr<argList> args;
        shared_ptr<Time> runTime;

        // Mesh
        autoPtr<dynamicFvMesh> meshPtr;
        dynamicFvMesh & mesh;

    protected:
        // Dictionaries
        IOdictionary couplingProperties;

        // FSI interface displacement
        wordList movingPatches;
        labelList movingPatchIDs;
        Field<vectorField> movingPatchesDispl;
        Field<vectorField> movingPatchesDisplOld;

        // Parallel FSI coupling variables
        labelList nPoints;
        labelList nGlobalPoints;
        labelList nGlobalCenters;
        labelList globalPointsUnique;
        labelList globalPointsNonUnique;

        float totalRunTime;
        int totalNbIterations;

        twoDPointCorrectorRBF twoDCorrector;

        std::unordered_map<unsigned int, std::unordered_map<std::string, unsigned int> > globalMovingPoints;
        std::vector<unsigned int> globalMovingPointLabels;

        labelList nbGlobalPoints;
};

#endif
