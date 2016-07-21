
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */


#ifndef foamSolidSolver_H
#define foamSolidSolver_H

#include <memory>

#include "BaseMultiLevelSolver.H"
#include "fvCFD.H"
#include "solidTractionFvPatchVectorField.H"
#include "leastSquaresVolPointInterpolation.H"
#include "twoDPointCorrectorRBF.H"
#include <time.h>
#include <deque>
#include <unordered_map>

using namespace fsi;
using std::shared_ptr;

class foamSolidSolver : public BaseMultiLevelSolver
{
    public:
        foamSolidSolver(
            string name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        ~foamSolidSolver();

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

        void getDisplacementLocal( matrix & displacement );

        void getReadPositionsLocal( matrix & readPositions );

        void getWritePositionsLocal( matrix & writePositions );

        void run();

        void setTractionLocal( const matrix & traction );

        string name;
        std::shared_ptr<argList> args;
        std::shared_ptr<Time> runTime;

        Foam::fvMesh mesh;

        static Foam::debug::debugSwitch debug;

    protected:
        void readCouplingProperties();

        // Fields
        volVectorField U;

        // Dictionaries
        IOdictionary couplingProperties;

        // Solid patch data
        wordList movingPatches;
        labelList movingPatchIDs;

        std::deque<solidTractionFvPatchVectorField *> forces;

        labelList nGlobalCenters;

        float totalRunTime;
        int totalNbIterations;

        twoDPointCorrectorRBF twoDCorrector;
        std::unordered_map<unsigned int, std::unordered_map<std::string, unsigned int> > movingPoints;
        std::vector<unsigned int> movingPointLabels;

        std::unordered_map<unsigned int, std::unordered_map<std::string, unsigned int> > globalMovingPoints;
        std::vector<unsigned int> globalMovingPointLabels;

        labelList nbGlobalPoints;
};

#endif
