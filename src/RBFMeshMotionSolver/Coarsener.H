
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <memory>
#include "ElRBFInterpolation.H"

namespace rbf {
class Coarsener {
    public:
        virtual ~Coarsener() {}

        virtual void compute(std::shared_ptr<RBFFunctionInterface> rbfFunction,
            std::unique_ptr<ElDistVector> positions,
            std::unique_ptr<ElDistVector> positionsInterpolation
            ) = 0;

        virtual bool initialized() = 0;

        virtual std::unique_ptr<ElDistVector> interpolate(const std::unique_ptr<ElDistVector> & values) = 0;
};
}  // namespace rbf
