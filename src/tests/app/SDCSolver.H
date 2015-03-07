
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCSolver_H
#define SDCSolver_H

#include <Eigen/Dense>

namespace sdc
{
  class SDCSolver
  {
public:

    virtual ~SDCSolver(){}

    virtual void evaluateFunction(
      const Eigen::VectorXd & q,
      const double t,
      Eigen::VectorXd & f
      ) = 0;

    virtual void implicitSolve(
      const double t,
      const double dt,
      const Eigen::VectorXd & qold,
      const Eigen::VectorXd rhs,
      Eigen::VectorXd & f,
      Eigen::VectorXd & result
      ) = 0;
  };
}

#endif
