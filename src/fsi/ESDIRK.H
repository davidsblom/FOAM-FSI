
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ESDIRK_H
#define ESDIRK_H

#include <memory>
#include "SDCSolver.H"
#include "fvCFD.H"

namespace sdc
{
    class ESDIRK
    {
public:

        ESDIRK(
            std::shared_ptr<SDCSolver> solver,
            std::string method
            );

        ~ESDIRK();

        void run();

        void solveTimeStep( const double t0 );

        std::shared_ptr<SDCSolver> solver;

private:

        bool isStageImplicit( double Akk );

        double dt;
        int nbStages;
        Eigen::MatrixXd A;
        Eigen::VectorXd B;
        Eigen::VectorXd C;
        int N;
    };
}

#endif
