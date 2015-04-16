
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
#include <time.h>

using namespace fsi;
using std::shared_ptr;

class foamFluidSolver : public BaseMultiLevelSolver
{
public:

    foamFluidSolver(
        string name,
        shared_ptr<argList> args,
        shared_ptr<Time> runTime
        );

    ~foamFluidSolver();

    virtual void getAcousticsDensityLocal( matrix & data ) = 0;

    virtual void getAcousticsVelocityLocal( matrix & data ) = 0;

    virtual void getAcousticsPressureLocal( matrix & data ) = 0;

    virtual void getTractionLocal( matrix & traction ) = 0;

    virtual void getWritePositionsLocalAcoustics( matrix & writePositions ) = 0;

    virtual void finalizeTimeStep();

    virtual void getReadPositions( matrix & readPositions );

    virtual void getWritePositions( matrix & writePositions );

    virtual void solve() = 0;

    virtual void solve(
        const matrix & input,
        matrix & output
        );

    int getInterfaceSize();

    int getInterfaceSizeLocal();

    void setDisplacementLocal( const matrix & displacement );

    void getReadPositionsLocal( matrix & readPositions );

    void getWritePositionsLocal( matrix & writePositions );

    void run();

    string name;
    shared_ptr<argList> args;
    shared_ptr<Time> runTime;

    // Mesh
    autoPtr<dynamicFvMesh> meshPtr;
    dynamicFvMesh & mesh;

protected:

    virtual void moveMesh( Field<vectorField> motion );

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
};

#endif
