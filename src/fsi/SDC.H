
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDC_H
#define SDC_H

#include <memory>
#include "SDCSolver.H"
#include "fvCFD.H"

namespace sdc
{
    class SDC
    {
public:

        SDC(
            std::shared_ptr<SDCSolver> solver,
            std::string rule,
            int nbNodes,
            double tol
            );

        ~SDC();

        void run();

        void solveTimeStep( const double t0 );

        std::shared_ptr<SDCSolver> solver;

        int nbNodes;

        Eigen::VectorXd nodes;
        Eigen::MatrixXd smat;
        Eigen::MatrixXd qmat;
        Eigen::VectorXd dsdc;
        double dt;
        int N;
        int k;
        double tol;
    };
}

#endif
