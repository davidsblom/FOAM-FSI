
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef IQNILSPostProcessing_H
#define IQNILSPostProcessing_H

#include <Eigen/QR>

#include "MultiLevelFsiSolver.H"
#include "PostProcessing.H"
#include "fvCFD.H"

namespace fsi
{
  class IQNILSPostProcessing : public PostProcessing
  {
public:

    IQNILSPostProcessing(
      shared_ptr<MultiLevelFsiSolver> fsi,
      int maxIter,
      double initialRelaxation,
      int maxUsedIterations,
      int nbReuse,
      double singularityLimit,
      int reuseInformationStartingFromTimeIndex,
      bool scaling
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

    virtual void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk,
      const matrix & B,
      const fsi::vector & xktilde,
      const fsi::vector & xkp
      );

    void removeColumnFromMatrix(
      matrix & A,
      int col
      );

    bool scaling;
    double singularityLimit;
    vector scalingFactors;
  };
}

#endif
