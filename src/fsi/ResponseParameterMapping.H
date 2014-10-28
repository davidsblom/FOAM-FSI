
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ResponseParameterMapping_H
#define ResponseParameterMapping_H

#include "SpaceMapping.H"

namespace fsi
{
  class ResponseParameterMapping : public SpaceMapping
  {
public:

    using SpaceMapping::performPostProcessing;

    ResponseParameterMapping(
      shared_ptr<SurrogateModel> fineModel,
      shared_ptr<SurrogateModel> coarseModel,
      int maxIter,
      int nbReuse,
      int reuseInformationStartingFromTimeIndex,
      double singularityLimit
      );

    virtual ~ResponseParameterMapping();

    virtual void iterationsConverged();

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

    std::deque<fsi::vector> xktildeList;
    std::deque<fsi::vector> xkList;
    std::deque<fsi::vector> mappedCoarseList;
    std::deque<std::deque<fsi::vector> > xktildePreviousTimeStepsList;
    std::deque<std::deque<fsi::vector> > xkPreviousTimeStepsList;
    std::deque<std::deque<fsi::vector> > mappedCoarsePreviousTimeStepsList;
  };
}

#endif
