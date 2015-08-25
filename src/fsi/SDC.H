
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDC_H
#define SDC_H

#include <memory>
#include "SDCSolver.H"
#include "fvCFD.H"
#include "AdaptiveTimeStepper.H"
#include "TimeIntegrationScheme.H"

namespace sdc
{
    class SDC : public TimeIntegrationScheme
    {
public:

        SDC(
            std::shared_ptr<SDCSolver> solver,
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
            std::string rule,
            int nbNodes,
            scalar tol,
            int minSweeps,
            int maxSweeps
            );

        SDC(
            std::string rule,
            int nbNodes,
            scalar tol
            );

        ~SDC();

        void solveTimeStep( const scalar t0 );

        virtual void run();

        virtual void getSourceTerm(
            const bool corrector,
            const int k,
            const scalar deltaT,
            fsi::vector & rhs,
            fsi::vector & qold
            );

        virtual void setFunction(
            const int k,
            const fsi::vector & f,
            const fsi::vector & result
            );

        virtual void setOldSolution(
            int timeIndex,
            const fsi::vector & result
            );

        virtual int getNbImplicitStages();

        virtual void outputResidual( std::string name );

        virtual bool isConverged();

        std::shared_ptr<SDCSolver> solver;
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper;

        int nbNodes;
        int N;
        int k;

        scalar dt;
        scalar tol;

        fsi::vector nodes;
        fsi::matrix smat;
        fsi::matrix qmat;
        fsi::vector nodesEmbedded;
        fsi::matrix smatEmbedded;
        fsi::matrix qmatEmbedded;
        fsi::vector dsdc;

private:

        void computeResidual(
            const fsi::matrix & qmat,
            const fsi::matrix & F,
            const scalar dt,
            fsi::matrix & qj
            );

        // Store function in memory in case the source term is requested
        // by the solver
        bool corrector;
        int stageIndex;
        fsi::matrix F;
        fsi::matrix Fold;
        fsi::matrix Sj;
        fsi::matrix solStages;
        bool convergence;
        int timeIndex;
        const int minSweeps;
        const int maxSweeps;
    };
}

#endif
