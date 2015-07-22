
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
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
            scalar singularityLimit
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

        void removeColumnFromMatrix(
            matrix & A,
            int col
            );
    };
}

#endif
