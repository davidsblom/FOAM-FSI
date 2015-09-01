
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
#include "SDCFsiSolverInterface.H"

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

using namespace fsi;
using std::shared_ptr;

class SDCFluidSolver : public foamFluidSolver, public sdc::SDCFsiSolverInterface
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
        const fsi::vector & q,
        const scalar t,
        fsi::vector & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual void getSolution( fsi::vector & solution, fsi::vector & f );

    virtual void setSolution(
        const fsi::vector & solution,
        const fsi::vector & f
        );

    virtual scalar getEndTime();

    virtual scalar getScalingFactor();

    virtual scalar getTimeStep();

    virtual scalar getStartTime();

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

    virtual void setDeltaT( scalar dt );

    virtual void setNumberOfImplicitStages( int k );

    virtual void nextTimeStep();

    virtual void getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        );

protected:

    void continuityErrs();

    void courantNo();

    void createFields();

    scalar evaluateMomentumResidual();

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
    volScalarField AU;
    volVectorField HU;
    volVectorField rhsU;
    surfaceScalarField rhsPhi;

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
    IOobject UFHeader;
    IOobject phiFHeader;
    volVectorField UF;
    surfaceScalarField phiF;

    bool turbulenceSwitch;
    bool explicitFirstStage;
};

#endif
