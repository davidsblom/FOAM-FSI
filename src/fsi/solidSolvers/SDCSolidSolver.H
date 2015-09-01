
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCSolidSolver_H
#define SDCSolidSolver_H

#include <memory>

#include "foamSolidSolver.H"
#include "constitutiveModel.H"
#include "RBFCoarsening.H"
#include "meshToMesh.H"
#include "SDCSolver.H"

using namespace fsi;
using namespace rbf;

class SDCSolidSolver : public foamSolidSolver, public sdc::SDCSolver
{
public:

    SDCSolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime
        );

    SDCSolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime,
        std::shared_ptr<rbf::RBFCoarsening> interpolator
        );

    virtual ~SDCSolidSolver();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void resetSolution();

    virtual void solve();

    virtual bool interpolateVolField( std::shared_ptr<BaseMultiLevelSolver> fineModel );

    virtual void evaluateFunction(
        const int k,
        const fsi::vector & q,
        const scalar t,
        fsi::vector & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual scalar getScalingFactor();

    virtual void getSolution( fsi::vector & solution );

    virtual void setSolution(
        const fsi::vector & solution,
        const fsi::vector & f
        );

    virtual scalar getEndTime();

    virtual scalar getTimeStep();

    virtual void nextTimeStep();

    virtual void setNumberOfImplicitStages( int k );

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

    virtual scalar getStartTime();

    virtual void getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        );

protected:

    void calculateEpsilonSigma();

    void initialize();

    void readSolidMechanicsControls();

    // Fields
    volTensorField gradU;
    volSymmTensorField epsilon;
    volSymmTensorField sigma;

    // Constitutive model
    constitutiveModel rheology;
    volScalarField rho;
    volScalarField mu;
    volScalarField lambda;
    surfaceScalarField muf;
    surfaceScalarField lambdaf;
    surfaceVectorField n;

    volVectorField V;

    // Solid mechanics controls
    const int minIter;
    const int maxIter;
    const scalar absoluteTolerance;
    const scalar relativeTolerance;

    const std::shared_ptr<rbf::RBFCoarsening> interpolator;

    int k;
    std::deque<volVectorField> UStages;
    std::deque<volVectorField> VStages;
    volVectorField rhsU;
    volVectorField rhsV;
    volVectorField UF;
    volVectorField VF;
};

#endif
