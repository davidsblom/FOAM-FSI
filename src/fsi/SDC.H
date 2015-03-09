
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDC_H
#define SDC_H

#include <memory>
#include "SDCSolver.H"

namespace sdc
{
  class SDC
  {
public:

    SDC(
      std::shared_ptr<SDCSolver> solver,
      int nbNodes
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
  };
}

#endif
