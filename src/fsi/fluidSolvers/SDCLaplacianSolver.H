
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
        const Eigen::VectorXd & q,
        const double t,
        Eigen::VectorXd & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual int getNbCells();

    virtual double getScalingFactor();

    virtual void getSolution( Eigen::VectorXd & solution );

    virtual void setSolution(
        const Eigen::VectorXd & solution,
        const Eigen::VectorXd & f
        );

    virtual double getEndTime();

    virtual double getTimeStep();

    virtual void nextTimeStep();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void setDeltaT( double dt );

    virtual void setNumberOfImplicitStages( int k );

    virtual void implicitSolve(
        bool corrector,
        const int k,
        const int kold,
        const double t,
        const double dt,
        const Eigen::VectorXd & qold,
        const Eigen::VectorXd & rhs,
        Eigen::VectorXd & f,
        Eigen::VectorXd & result
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
