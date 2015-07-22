
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *
 * Adaptive time stepping based on :
 * Soderlind and Wang, Adaptive time-stepping and computational stability,
 * Journal of Computational and Applied Mathematics 185 (2006) 225-243.
 */

#ifndef AdaptiveTimeStepper_H
#define AdaptiveTimeStepper_H

#include <memory>
#include <assert.h>
#include <math.h>
#include "BaseMultiLevelSolver.H"
#include "fvCFD.H"

namespace sdc
{
    class AdaptiveTimeStepper
    {
public:

        AdaptiveTimeStepper( bool enabled );

        AdaptiveTimeStepper(
            bool enabled,
            std::string filter,
            scalar tol,
            scalar safetyFactor
            );

        ~AdaptiveTimeStepper();

        bool determineNewTimeStep(
            const fsi::vector & errorEstimate,
            const fsi::vector & solution,
            const scalar computedTimeStep,
            scalar & newTimeStep
            );

        bool determineNewTimeStep(
            const scalar errorEstimate,
            const scalar computedTimeStep,
            scalar & newTimeStep
            );

        scalar elementary(
            const scalar c1,
            const scalar c0,
            const scalar rho
            );

        scalar h211b(
            const scalar c1,
            const scalar c0,
            const scalar rho
            );

        bool isAccepted();

        bool isEnabled();

        bool isPreviousStepAccepted();

        scalar limit( const scalar u );

        scalar pi42(
            const scalar c1,
            const scalar c0,
            const scalar rho
            );

        void setEndTime( scalar endTime );

        void setOrderEmbeddedMethod( int order );

private:

        bool enabled;
        std::string filter;
        scalar tol;
        scalar safetyFactor;
        int k;
        scalar cerrold;
        scalar rhoold;
        int timeStepIndex;
        bool accepted;
        bool previousTimeStepRejected;
        scalar endTime;
        scalar t;
    };
}

#endif
