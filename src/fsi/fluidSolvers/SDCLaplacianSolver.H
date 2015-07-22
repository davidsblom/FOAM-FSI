
/*
 * Class
 *   SDCLaplacianSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCLaplacianSolver_H
#define SDCLaplacianSolver_H

#include <memory>
#include <deque>
#include "fvCFD.H"
#include "SDCSolver.H"

using std::shared_ptr;

class SDCLaplacianSolver : public sdc::SDCSolver
{
public:

    SDCLaplacianSolver(
        string name,
        shared_ptr<argList> args,
        shared_ptr<Time> runTime
        );

    virtual ~SDCLaplacianSolver();

    virtual void evaluateFunction(
        const int k,
        const fsi::vector & q,
        const scalar t,
        fsi::vector & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual int getNbCells();

    virtual scalar getScalingFactor();

    virtual void getSolution( fsi::vector & solution );

    virtual void setSolution(
        const fsi::vector & solution,
        const fsi::vector & f
        );

    virtual scalar getEndTime();

    virtual scalar getTimeStep();

    virtual void nextTimeStep();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void setDeltaT( scalar dt );

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

private:

    bool init;

    string name;
    std::shared_ptr<argList> args;
    std::shared_ptr<Time> runTime;

    Foam::fvMesh mesh;

    volScalarField T;
    volScalarField rhsT;
    IOdictionary transportProperties;
    dimensionedScalar DT;
    int nNonOrthCorr;

    int k;
    std::deque<volScalarField> TStages;
};

#endif
