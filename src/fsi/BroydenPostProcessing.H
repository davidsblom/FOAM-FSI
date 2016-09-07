
/*
 * Copyright [2016] <David Blom>
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
                int maxIter,
                scalar initialRelaxation,
                int maxUsedIterations,
                int nbReuse,
                scalar singularityLimit,
                int reuseInformationStartingFromTimeIndex
                );

            virtual ~BroydenPostProcessing();

            virtual void finalizeStage();

            virtual void finalizeTimeStep();

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

            matrix J;
            const scalar singularityLimit;
    };
} // namespace fsi

#endif
