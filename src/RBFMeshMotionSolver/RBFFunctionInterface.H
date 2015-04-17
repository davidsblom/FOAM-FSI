
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFFunctionInterface_H
#define RBFFunctionInterface_H

#include <math.h>
#include <string>
#include <memory>
#include <Eigen/Dense>

namespace rbf
{
    class RBFFunctionInterface
    {
public:

        virtual ~RBFFunctionInterface(){}

        virtual double evaluate( double value ) = 0;
    };
}

#endif
