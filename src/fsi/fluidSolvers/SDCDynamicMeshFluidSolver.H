
/*
 * Class
 *   SDCDynamicMeshFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCDynamicMeshFluidSolver_H
#define SDCDynamicMeshFluidSolver_H

#include "SDCMovingWallVelocityFvPatchVectorField.H"
#include "SDCFluidSolver.H"
#include "SDCFsiSolverInterface.H"

using namespace fsi;
using std::shared_ptr;

class SDCDynamicMeshFluidSolver : public SDCFluidSolver
{
    public:
        SDCDynamicMeshFluidSolver(
            std::string name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        virtual ~SDCDynamicMeshFluidSolver();

        virtual void resetSolution();

        virtual void evaluateFunction(
            const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            );

        virtual int getDOF();

        virtual void getSolution(
            fsi::vector & solution,
            fsi::vector & f
            );

        virtual void setSolution(
            const fsi::vector & solution,
            const fsi::vector & f
            );

        virtual void implicitSolve(
            bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs,
            fsi::vector & f,
            fsi::vector & result
            );

        virtual void prepareImplicitSolve(
            bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs
            );

        virtual void setNumberOfImplicitStages( int k );

        virtual void nextTimeStep();

        virtual void getVariablesInfo(
            std::deque<int> & dof,
            std::deque<bool> & enabled,
            std::deque<std::string> & names
            );

        virtual bool isConverged();

        virtual void solve();

        virtual void getTractionLocal( matrix & traction );

    private:
        void createFields();

        scalar evaluateMomentumResidual();

        void initialize();

        // Fields
        surfaceVectorField Uf;
        surfaceVectorField rhsUf;
        surfaceScalarField rhsMeshPhi;

        std::deque<surfaceVectorField> UfStages;
        std::deque<pointField> pointsStages;
        std::deque<volScalarField::DimensionedInternalField> volumeStages;
        std::deque<surfaceScalarField> interpolateVolumeStages;
        IOobject UfFHeader;
        surfaceVectorField UfF;
        surfaceScalarField meshPhiF;

        int kold;
        int indexk;
        int sweep;
};

#endif
