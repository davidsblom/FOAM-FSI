
/*
 * Class
 *   SDCFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCFluidSolver_H
#define SDCFluidSolver_H

#include <memory>
#include <deque>

#include "foamFluidSolver.H"
#include "SDCSolver.H"

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "faceTetPolyPatch.H"
#include "fixedValuePointPatchFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "tetPolyPatchInterpolation.H"

using namespace fsi;
using std::shared_ptr;

class SDCFluidSolver : public foamFluidSolver, public sdc::SDCSolver
{
public:

    SDCFluidSolver(
        string name,
        shared_ptr<argList> args,
        shared_ptr<Time> runTime
        );

    virtual ~SDCFluidSolver();

    virtual void getAcousticsDensityLocal( matrix & data );

    virtual void getAcousticsVelocityLocal( matrix & data );

    virtual void getAcousticsPressureLocal( matrix & data );

    virtual void getTractionLocal( matrix & traction );

    virtual void getWritePositionsLocalAcoustics( matrix & writePositions );

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void resetSolution();

    virtual void solve();

    virtual void evaluateFunction(
        const int k,
        const Eigen::VectorXd & q,
        const double t,
        Eigen::VectorXd & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual void getSolution( Eigen::VectorXd & solution );

    virtual int getNbCells();

    virtual double getScalingFactor();

    virtual double getTimeStep();

    virtual void implicitSolve(
        bool corrector,
        const int k,
        const double t,
        const double dt,
        const Eigen::VectorXd & qold,
        const Eigen::VectorXd & rhs,
        Eigen::VectorXd & f,
        Eigen::VectorXd & result
        );

    virtual void setDeltaT( double dt );

    virtual void setNumberOfStages( int k );

    virtual void nextTimeStep();

private:

    void continuityErrs();

    void courantNo();

    void createFields();

    double evaluateMomentumResidual( double dt );

    void initialize();

    void readPIMPLEControls();

    // Dictionaries
    IOdictionary transportProperties;
    dictionary pimple;

    // Dimensioned scalars
    dimensionedScalar nu;
    dimensionedScalar rho;

    // Fields
    volScalarField p;
    volVectorField U;
    surfaceScalarField phi;
    surfaceVectorField Uf;
    volScalarField AU;
    volVectorField HU;
    volVectorField rhsU;
    surfaceVectorField rhsUf;

    // PIMPLE controls
    int nCorr;
    int nNonOrthCorr;
    int minIter;
    int maxIter;
    scalar absoluteTolerance;
    scalar relativeTolerance;

    // Continuity errors
    scalar sumLocalContErr;
    scalar globalContErr;
    scalar cumulativeContErr;

    label pRefCell;
    scalar pRefValue;
    singlePhaseTransportModel laminarTransport;
    autoPtr<incompressible::turbulenceModel> turbulence;

    int k;
    std::deque<volScalarField> pStages;
    std::deque<surfaceScalarField> phiStages;
    std::deque<volVectorField> UStages;
    std::deque<surfaceVectorField> UfStages;
    volVectorField UF;
    surfaceVectorField UfF;

    bool turbulenceSwitch;
};

#endif
