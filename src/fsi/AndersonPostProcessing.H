
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef AndersonPostProcessing_H
#define AndersonPostProcessing_H

#include <Eigen/QR>

#include "MultiLevelFsiSolver.H"
#include "PostProcessing.H"
#include "fvCFD.H"

namespace fsi
{
  class AndersonPostProcessing : public PostProcessing
  {
public:

    AndersonPostProcessing(
      shared_ptr<MultiLevelFsiSolver> fsi,
      int maxIter,
      double initialRelaxation,
      int maxUsedIterations,
      int nbReuse,
      double singularityLimit,
      int reuseInformationStartingFromTimeIndex,
      bool scaling,
      double beta,
      bool updateJacobian
      );

    void fixedUnderRelaxation(
      vector & xk,
      vector & R,
      vector & yk
      );

    virtual void performPostProcessing(
      const vector & x0,
      vector & xk
      );

    virtual void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk
      );

    virtual void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk,
      bool residualCriterium
      );

    void applyScaling( vector & vec );

    void applyScaling( matrix & mat );

    void removeColumnFromMatrix(
      matrix & A,
      int col
      );

    void determineScalingFactors( const vector & output );

    bool scaling;
    double beta;
    double singularityLimit;
    bool updateJacobian;
    vector scalingFactors;
    matrix Jprev;
    int sizeVar0;
    int sizeVar1;
  };
}

#endif
