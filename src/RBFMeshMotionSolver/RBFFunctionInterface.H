
/*
 * Copyright [2016] <David Blom>
 */

#ifndef RBFFunctionInterface_H
#define RBFFunctionInterface_H

#include <math.h>
#include <Eigen/Dense>
#include <string>
#include <memory>
#include "fvCFD.H"

namespace rbf {
class RBFFunctionInterface {
    public:
        virtual ~RBFFunctionInterface() {}

        virtual scalar evaluate(scalar value) = 0;
};
}  // namespace rbf

#endif
