
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
#include "fvCFD.H"

namespace sdc
{
    class AdaptiveTimeStepper
    {
public:

        AdaptiveTimeStepper(
            bool enabled,
            std::string filter,
            double tol,
            int k
            );

        ~AdaptiveTimeStepper();

        bool determineNewTimeStep(
            const double errorEstimate,
            const double computedTimeStep,
            double & newTimeStep
            );

        double elementary(
            const double c1,
            const double c0,
            const double rho
            );

        double h211b(
            const double c1,
            const double c0,
            const double rho
            );

        double limit( const double u );

        double pi42(
            const double c1,
            const double c0,
            const double rho
            );

        bool enabled;
        std::string filter;
        double tol;
        int k;
        double cerrold;
        double rhoold;
        int timeStepIndex;
        bool accepted;
    };
}

#endif
