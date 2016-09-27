
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
            const std::string & name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        virtual ~SDCFluidSolver();

        virtual void getAcousticsVelocityLocal( matrix & data ) override;

        virtual void getAcousticsPressureLocal( matrix & data ) override;

        virtual void getAcousticsDensityLocal( matrix & data ) override;

        virtual void getTractionLocal( matrix & traction ) override;

        virtual void getWritePositionsLocalAcoustics( matrix & writePositions ) override;

        virtual void initTimeStep() override;

        virtual bool isRunning() override;

        virtual void resetSolution() override;

        virtual void setAcousticsTemperatureGradient( const matrix & temperatureGradient ) override;

        virtual void setAcousticsPressureGradient( const matrix & pressureGradient ) override;

        virtual void setAcousticsVelocityGradient( const matrix & velocityGradient ) override;

        virtual void solve() override;

        virtual void evaluateFunction(
            const int k,
            const fsi::vector & q,
            const scalar t,
            fsi::vector & f
            ) override;

        virtual void finalizeTimeStep() override;

        virtual int getDOF() override;

        virtual void getSolution(
            fsi::vector & solution,
            fsi::vector & f
            ) override;

        virtual void setSolution(
            const fsi::vector & solution,
            const fsi::vector & f
            ) override;

        virtual scalar getEndTime() override;

        virtual scalar getTimeStep() override;

        virtual scalar getStartTime() override;

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
            ) override;

        void prepareImplicitSolve(
            bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs
            ) override;

        virtual void setDeltaT( scalar dt );

        virtual void setNumberOfImplicitStages( int k ) override;

        virtual void nextTimeStep() override;

        virtual void getVariablesInfo(
            std::deque<int> & dof,
            std::deque<bool> & enabled,
            std::deque<std::string> & names
            ) override;

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
        label writeInterval;
};

#endif
