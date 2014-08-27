
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef BroydenPostProcessing_H
#define BroydenPostProcessing_H

#include "PostProcessing.H"

namespace fsi
{
  class BroydenPostProcessing : public PostProcessing
  {
public:

    BroydenPostProcessing(
      shared_ptr<MultiLevelFsiSolver> fsi,
      double initialRelaxation,
      int maxIter,
      int maxUsedIterations,
      int nbReuse,
      int reuseInformationStartingFromTimeIndex
      );

    virtual ~BroydenPostProcessing();

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
  };
}

#endif
