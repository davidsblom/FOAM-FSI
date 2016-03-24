
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "BaseMultiLevelSolver.H"
#include <deque>

namespace sdc
{
    class SDCSolver
    {
        public:
            virtual ~SDCSolver(){}

            virtual void evaluateFunction(
                const int k,
                const fsi::vector & q,
                const scalar t,
                fsi::vector & f
                ) = 0;

            virtual void finalizeTimeStep() = 0;

            virtual int getDOF() = 0;

            virtual void getSolution(
                fsi::vector & solution,
                fsi::vector & f
                ) = 0;

            virtual void setSolution(
                const fsi::vector & solution,
                const fsi::vector & f
                ) = 0;

            virtual scalar getEndTime() = 0;

            virtual scalar getTimeStep() = 0;

            virtual void nextTimeStep() = 0;

            virtual void initTimeStep() = 0;

            virtual void setNumberOfImplicitStages( int k ) = 0;

            virtual void implicitSolve(
                bool corrector,
                const int k,
                const int kold,
                const scalar t,
                const scalar dt,
                const fsi::vector & qold,
                const fsi::vector & rhs,
                fsi::vector & f,
                fsi::vector & result
                ) = 0;

            virtual scalar getStartTime()
            {
                return 0;
            }

            virtual void getVariablesInfo(
                std::deque<int> & dof,
                std::deque<bool> & enabled,
                std::deque<std::string> & names
                )
            {
                dof.push_back( 1 );
                enabled.push_back( true );
                names.push_back( "0" );
            }

            virtual bool isConverged()
            {
                return true;
            }
    };
}
