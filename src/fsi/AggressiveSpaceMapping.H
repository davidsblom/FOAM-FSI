
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef AggressiveSpaceMapping_H
#define AggressiveSpaceMapping_H

#include "SpaceMapping.H"

namespace fsi
{
  class AggressiveSpaceMapping : public SpaceMapping
  {
public:

    using SpaceMapping::performPostProcessing;

    AggressiveSpaceMapping(
      shared_ptr<SurrogateModel> fineModel,
      shared_ptr<SurrogateModel> coarseModel,
      int maxIter,
      int maxUsedIterations,
      int nbReuse,
      int reuseInformationStartingFromTimeIndex,
      double singularityLimit
      );

    virtual ~AggressiveSpaceMapping();

    virtual void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk,
      bool residualCriterium
      );
  };
}

#endif
