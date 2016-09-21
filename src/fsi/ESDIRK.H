
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef ESDIRK_H
#define ESDIRK_H

#include <memory>
#include "SDCSolver.H"
#include "AdaptiveTimeStepper.H"
#include "fvCFD.H"
#include "TimeIntegrationScheme.H"

namespace sdc
{
    class ESDIRK : public TimeIntegrationScheme
    {
        public:
            ESDIRK(
                std::shared_ptr<SDCSolver> solver,
                const std::string & method,
                std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper
                );

            explicit ESDIRK( const std::string & method );

            ~ESDIRK();

            void run();

            void solveTimeStep( const scalar t0 );

            virtual void getSourceTerm(
                const bool corrector,
                const int k,
                const int sweep,
                const scalar deltaT,
                fsi::vector & rhs,
                fsi::vector & qold
                );

            virtual void setFunction(
                const int k,
                const fsi::vector & f,
                const fsi::vector & result
                );

            virtual int getNbImplicitStages();

            virtual void setOldSolution(
                int timeIndex,
                const fsi::vector & result
                );

            virtual void outputResidual( const std::string & name );

            virtual bool isConverged();

            static bool isStageImplicit( scalar Akk );

            std::shared_ptr<SDCSolver> solver;
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper;
            scalar dt;
            int nbStages;
            fsi::matrix A;
            fsi::vector B;
            fsi::vector C;
            fsi::vector Bhat;
            int N;

        private:
            void initializeButcherTableau( const std::string & method );

            // Store function in memory in case the source term is requested
            // by the solver. This functionality is used by the moving wall
            // velocity boundary condition.
            int stageIndex;
            fsi::matrix F;
            fsi::matrix solStages;
            fsi::vector qold;
            int timeIndex;
    };
}

#endif
