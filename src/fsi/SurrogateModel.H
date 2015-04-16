
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SurrogateModel_H
#define SurrogateModel_H

#include "MultiLevelSolver.H"

namespace fsi
{
    class SurrogateModel
    {
public:

        virtual ~SurrogateModel(){}

        virtual bool allConverged() = 0;

        virtual void evaluate(
            const fsi::vector & x,
            fsi::vector & output,
            fsi::vector & R
            ) = 0;

        virtual bool isConvergence() = 0;

        virtual bool isConvergence(
            const fsi::vector & x,
            const fsi::vector & xprev
            ) = 0;

        virtual void optimize(
            const fsi::vector & x0,
            fsi::vector & xk
            ) = 0;

        virtual void optimize(
            const fsi::vector & y,
            const fsi::vector & x0,
            fsi::vector & xk
            ) = 0;

        virtual void setSurrogateData(
            fsi::vector & xf,
            matrix & J
            ) = 0;

        virtual void setUseJacobian( bool useJacobian ) = 0;
    };
}

#endif
