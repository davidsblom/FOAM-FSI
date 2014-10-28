
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ManifoldMapping_H
#define ManifoldMapping_H

#include "SpaceMapping.H"

namespace fsi
{
  class ManifoldMapping : public SpaceMapping
  {
public:

    using SpaceMapping::performPostProcessing;

    ManifoldMapping(
      shared_ptr<SurrogateModel> fineModel,
      shared_ptr<SurrogateModel> coarseModel,
      int maxIter,
      int nbReuse,
      int reuseInformationStartingFromTimeIndex,
      double singularityLimit
      );

    virtual ~ManifoldMapping();

    virtual void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk,
      bool residualCriterium
      );

    void removeColumnFromMatrix(
      matrix & A,
      int col
      );

    int iter;
  };
}

#endif
