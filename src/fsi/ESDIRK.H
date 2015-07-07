
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ESDIRK_H
#define ESDIRK_H

#include <memory>
#include "SDCSolver.H"
#include "AdaptiveTimeStepper.H"
#include "fvCFD.H"

namespace sdc
{
    class ESDIRK
    {
public:

        ESDIRK(
            std::shared_ptr<SDCSolver> solver,
            std::string method,
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper
            );

        ~ESDIRK();

        void run();

        void solveTimeStep( const double t0 );

private:

        bool isStageImplicit( double Akk );

        std::shared_ptr<SDCSolver> solver;
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper;
        double dt;
        int nbStages;
        Eigen::MatrixXd A;
        Eigen::VectorXd B;
        Eigen::VectorXd C;
        Eigen::VectorXd Bhat;
        int N;
    };
}

#endif
