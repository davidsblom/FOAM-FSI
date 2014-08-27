
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef OutputSpaceMapping_H
#define OutputSpaceMapping_H

#include <deque>

#include "SurrogateModel.H"
#include "SpaceMapping.H"

using std::deque;

namespace fsi
{
  class OutputSpaceMapping : public SpaceMapping
  {
public:

    using SpaceMapping::performPostProcessing;

    OutputSpaceMapping(
      shared_ptr<SurrogateModel> fineModel,
      shared_ptr<SurrogateModel> coarseModel,
      int maxIter,
      int nbReuse,
      int reuseInformationStartingFromTimeIndex
      );

    virtual ~OutputSpaceMapping();

    virtual void iterationsConverged();

    virtual void performPostProcessing(
      const vector & y,
      const vector & x0,
      vector & xk,
      bool residualCriterium
      );

    shared_ptr<SurrogateModel> surrogateModel;
    deque<vector> sols;
    deque<deque<vector> > solsList;
  };
}

#endif
