
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "SDCSolver.H"

namespace sdc
{
    class SDCFsiSolverInterface : public SDCSolver
    {
public:

        virtual ~SDCFsiSolverInterface(){}

        virtual void prepareImplicitSolve(
            bool corrector,
            const int k,
            const int kold,
            const scalar t,
            const scalar dt,
            const fsi::vector & qold,
            const fsi::vector & rhs
            ) = 0;

        virtual void finalizeImplicitSolve( int /*k*/ ){}
    };
}
