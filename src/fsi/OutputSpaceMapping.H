
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
            int maxUsedIterations,
            int nbReuse,
            int reuseInformationStartingFromTimeIndex,
            scalar singularityLimit,
            int order
            );

        virtual ~OutputSpaceMapping();

        virtual void finalizeTimeStep();

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

        shared_ptr<SurrogateModel> surrogateModel;
        deque<vector> sols;
        deque<deque<vector> > solsList;
        deque<deque<deque<vector> > > solsTimeList;

        const int order;
    };
}

#endif
