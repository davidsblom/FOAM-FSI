
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SpaceMapping_H
#define SpaceMapping_H

#include <deque>

#include "SurrogateModel.H"
#include "fvCFD.H"

using namespace fsi;
using std::deque;

class SpaceMapping
{
    public:
        SpaceMapping(
            shared_ptr<SurrogateModel> fineModel,
            shared_ptr<SurrogateModel> coarseModel,
            int maxIter,
            int maxUsedIterations,
            int nbReuse,
            int reuseInformationStartingFromTimeIndex,
            scalar singularityLimit
            );

        virtual ~SpaceMapping(){}

        virtual void finalizeTimeStep();

        virtual bool isConvergence(
            const fsi::vector & xk,
            const fsi::vector & xkprev,
            bool residualCriterium
            );

        virtual void iterationsConverged();

        virtual void performPostProcessing(
            const fsi::vector & x0,
            fsi::vector & xk
            );

        virtual void performPostProcessing(
            const fsi::vector & y,
            const fsi::vector & x0,
            fsi::vector & xk
            );

        virtual void performPostProcessing(
            const fsi::vector & y,
            const fsi::vector & x0,
            fsi::vector & xk,
            bool residualCriterium
            ) = 0;

        shared_ptr<SurrogateModel> fineModel;
        shared_ptr<SurrogateModel> coarseModel;

        const int maxIter;
        const int maxUsedIterations;
        const int nbReuse;
        const int reuseInformationStartingFromTimeIndex;
        const scalar singularityLimit;
        int timeIndex;

        deque<fsi::vector> coarseResiduals;
        deque<fsi::vector> fineResiduals;
        deque<deque<fsi::vector> > coarseResidualsList;
        deque<deque<fsi::vector> > fineResidualsList;
        deque<deque<deque<fsi::vector> > > coarseResidualsTimeList;
        deque<deque<deque<fsi::vector> > > fineResidualsTimeList;
};

#endif
