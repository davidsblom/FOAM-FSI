
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef AitkenPostProcessing_H
#define AitkenPostProcessing_H

#include "PostProcessing.H"

namespace fsi
{
  class AitkenPostProcessing : public PostProcessing
  {
public:

    AitkenPostProcessing(
      shared_ptr<MultiLevelFsiSolver> fsi,
      double initialRelaxation,
      int maxIter,
      int maxUsedIterations,
      int nbReuse,
      int reuseInformationStartingFromTimeIndex
      );

    virtual ~AitkenPostProcessing();

    void performPostProcessing(
      const vector & x0,
      vector & xk
      );

    void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk
      );

    void performPostProcessing(
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

    double aitkenFactor;
  };
}

#endif
