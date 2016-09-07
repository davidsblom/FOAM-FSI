
/*
 * Copyright [2016] <David Blom>
 */

#ifndef ASMILS_H
#define ASMILS_H

#include "SpaceMapping.H"

namespace fsi
{
    class ASMILS : public SpaceMapping
    {
        public:
            using SpaceMapping::performPostProcessing;

            ASMILS(
                shared_ptr<SurrogateModel> fineModel,
                shared_ptr<SurrogateModel> coarseModel,
                int maxIter,
                int maxUsedIterations,
                int nbReuse,
                int reuseInformationStartingFromTimeIndex,
                scalar singularityLimit,
                scalar beta
                );

            virtual ~ASMILS();

            virtual void performPostProcessing(
                const vector & y,
                const vector & x0,
                vector & xk,
                bool residualCriterium
                );

            void fixedUnderRelaxation(
                vector & xk,
                vector & zstar,
                vector & zk
                );

        private:
            scalar beta;
    };
} // namespace fsi

#endif
