
/*
 * Copyright [2016] <David Blom>
 */

#ifndef AitkenPostProcessing_H
#define AitkenPostProcessing_H

#include <memory>
#include "PostProcessing.H"

namespace fsi {
class AitkenPostProcessing : public PostProcessing {
    public:
        AitkenPostProcessing(shared_ptr<MultiLevelFsiSolver> fsi,
            scalar initialRelaxation,
            int maxIter,
            int maxUsedIterations,
            int nbReuse,
            int reuseInformationStartingFromTimeIndex
            );

        virtual ~AitkenPostProcessing();

        void performPostProcessing(const vector & x0,
            vector & xk
            );

        void performPostProcessing(const vector & y,
            const vector & x0,
            vector & xk
            );

        void performPostProcessing(const vector & y,
            const vector & x0,
            vector & xk,
            bool residualCriterium
            );

        scalar aitkenFactor;
};
} // namespace fsi

#endif
