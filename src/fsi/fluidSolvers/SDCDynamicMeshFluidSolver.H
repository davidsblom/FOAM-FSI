
/*
 * Class
 *   SDCDynamicMeshFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCDynamicMeshFluidSolver_H
#define SDCDynamicMeshFluidSolver_H

#include <memory>
#include <deque>

#include "foamFluidSolver.H"
#include "SDCSolver.H"

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "SDCMovingWallVelocityFvPatchVectorField.H"

using namespace fsi;
using std::shared_ptr;

class SDCDynamicMeshFluidSolver : public foamFluidSolver, public sdc::SDCSolver
{
public:

    SDCDynamicMeshFluidSolver(
        string name,
        shared_ptr<argList> args,
        shared_ptr<Time> runTime
        );

    virtual ~SDCDynamicMeshFluidSolver();

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

    virtual void setSolution(
        const Eigen::VectorXd & solution,
        const Eigen::VectorXd & f
        );

    virtual double getEndTime();

    virtual double getScalingFactor();

    virtual double getTimeStep();

    virtual double getStartTime();

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

    virtual void setDeltaT( double dt );

    virtual void setNumberOfStages( int k );

    virtual void nextTimeStep();

    virtual void getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        );

    virtual bool isConverged();

private:

    void continuityErrs();

    void createFields();

    double evaluateMomentumResidual();

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
    surfaceVectorField Uf;
    surfaceScalarField phi;
    volScalarField AU;
    volVectorField HU;
    volVectorField rhsU;
    surfaceVectorField rhsUf;
    surfaceScalarField rhsMeshPhi;

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
    std::deque<pointField> pointsStages;
    std::deque<volScalarField::DimensionedInternalField> volumeStages;
    std::deque<surfaceScalarField> interpolateVolumeStages;
    IOobject UFHeader;
    IOobject phiFHeader;
    volVectorField UF;
    surfaceVectorField UfF;
    surfaceScalarField meshPhiF;

    bool turbulenceSwitch;
    bool explicitFirstStage;
};

#endif
