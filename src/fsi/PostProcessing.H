
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef PostProcessing_H
#define PostProcessing_H

#include "MultiLevelFsiSolver.H"

using std::deque;

namespace fsi
{
    class PostProcessing
    {
public:

        PostProcessing(
            shared_ptr<MultiLevelFsiSolver> fsi,
            scalar initialRelaxation,
            int maxIter,
            int maxUsedIterations,
            int nbReuse,
            int reuseInformationStartingFromTimeIndex
            );

        virtual ~PostProcessing(){}

        virtual void performPostProcessing(
            const vector & x0,
            vector & xk
            ) = 0;

        virtual void performPostProcessing(
            const vector & y,
            const vector & x0,
            vector & xk
            ) = 0;

        virtual void performPostProcessing(
            const vector & y,
            const vector & x0,
            vector & xk,
            bool residualCriterium
            ) = 0;

        virtual void finalizeTimeStep();

        bool isConvergence(
            const vector & xk,
            const vector & xkprev,
            bool residualCriterium
            );

        void iterationsConverged( bool keepIterations );

        shared_ptr<MultiLevelFsiSolver> fsi;

        scalar initialRelaxation;

        const int maxIter;
        const int maxUsedIterations;
        const int nbReuse;
        int timeIndex;
        const int reuseInformationStartingFromTimeIndex;

        deque<vector> residuals;
        deque<vector> sols;
        deque<deque<vector> > residualsList;
        deque<deque<vector> > solsList;
        deque<deque<deque<vector> > > residualsTimeList;
        deque<deque<deque<vector> > > solsTimeList;
    };
}

#endif
