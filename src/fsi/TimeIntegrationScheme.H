
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TimeIntegrationScheme_H
#define TimeIntegrationScheme_H

#include "DataValues.H"

namespace sdc
{
    class TimeIntegrationScheme
    {
        public:
            virtual ~TimeIntegrationScheme(){}

            virtual void getSourceTerm(
                const bool corrector,
                const int k,
                const int sweep,
                const scalar deltaT,
                fsi::vector & rhs,
                fsi::vector & qold
                ) = 0;

            virtual void setFunction(
                const int k,
                const fsi::vector & f,
                const fsi::vector & result
                ) = 0;

            virtual int getNbImplicitStages() = 0;

            virtual void setOldSolution(
                int timeIndex,
                const fsi::vector & result
                ) = 0;

            virtual void outputResidual( std::string name ) = 0;

            virtual bool isConverged() = 0;

            virtual void run() = 0;
    };
}

#endif
