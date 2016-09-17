
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <memory>
#include "Coarsener.H"

namespace rbf {
class NoCoarsening : public Coarsener {
    public:
        NoCoarsening();

        ~NoCoarsening();

        void compute(std::shared_ptr<RBFFunctionInterface> rbfFunction,
            std::unique_ptr<ElDistVector> positions,
            std::unique_ptr<ElDistVector> positionsInterpolation
            );

        bool initialized();

        std::unique_ptr<ElDistVector> interpolate(const std::unique_ptr<ElDistVector> & values);

    private:
        std::unique_ptr<ElRBFInterpolation> rbf;
};
}  // namespace rbf
