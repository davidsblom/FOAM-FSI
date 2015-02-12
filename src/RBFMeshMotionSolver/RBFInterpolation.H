
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include <string>
#include <memory>
#include <Eigen/Dense>
#include "RBFFunctionInterface.H"

namespace rbf
{
  typedef Eigen::MatrixXd matrix;
  typedef Eigen::VectorXd vector;

  class RBFInterpolation
  {
public:

    RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction );

    void compute(
      const matrix & positions,
      const matrix & positionsInterpolation
      );

    void interpolate(
      const matrix & values,
      matrix & valuesInterpolation
      );

    void interpolate(
      const matrix & positions,
      const matrix & positionsInterpolation,
      const matrix & values,
      matrix & valuesInterpolation
      );

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    bool computed;
    int n_A;
    int n_B;
    int dimGrid;
    matrix Hhat;
    matrix Phi;
  };
}

#endif
