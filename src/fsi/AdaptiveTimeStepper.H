
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
#include <Eigen/Dense>

namespace sdc
{
    class AdaptiveTimeStepper
    {
public:

        AdaptiveTimeStepper(
            bool enabled,
            std::string filter,
            double tol,
            double safetyFactor
            );

        ~AdaptiveTimeStepper();

        bool determineNewTimeStep(
            const Eigen::VectorXd & errorEstimate,
            const Eigen::VectorXd & solution,
            const double computedTimeStep,
            double & newTimeStep
            );

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

        bool isAccepted();

        bool isEnabled();

        bool isPreviousStepAccepted();

        double limit( const double u );

        double pi42(
            const double c1,
            const double c0,
            const double rho
            );

        void setEndTime( double endTime );

        void setOrderEmbeddedMethod( int order );

private:

        bool enabled;
        std::string filter;
        double tol;
        double safetyFactor;
        int k;
        double cerrold;
        double rhoold;
        int timeStepIndex;
        bool accepted;
        bool previousTimeStepRejected;
        double endTime;
        double t;
    };
}

#endif
