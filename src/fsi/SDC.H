
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDC_H
#define SDC_H

#include <memory>
#include "SDCSolver.H"
#include "fvCFD.H"
#include "TimeIntegrationScheme.H"
#include "QuadratureInterface.H"
#include "DataStorage.H"

namespace sdc
{
    class SDC : public TimeIntegrationScheme
    {
        public:
            SDC(
                std::shared_ptr<SDCSolver> solver,
                std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
                std::shared_ptr<sdc::DataStorage> data,
                scalar tol,
                int minSweeps,
                int maxSweeps
                );

            SDC(
                std::shared_ptr<SDCSolver> solver,
                std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
                scalar tol,
                int minSweeps,
                int maxSweeps
                );

            SDC(
                std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
                scalar tol
                );

            ~SDC();

            void init();

            void solveTimeStep( const scalar t0 );

            virtual void run();

            void computeResidual(
                const fsi::matrix & qmat,
                const scalar dt,
                fsi::matrix & qj
                );

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

            virtual void setOldSolution(
                int timeIndex,
                const fsi::vector & result
                );

            virtual int getNbImplicitStages();

            virtual void outputResidual( const std::string & name );

            virtual bool isConverged();

            std::shared_ptr<SDCSolver> solver;

            int nbNodes;
            int N;
            int k;
            int sweep;

            scalar dt;
            scalar tol;

            fsi::vector nodes;
            fsi::matrix smat;
            fsi::matrix qmat;
            fsi::vector nodesEmbedded;
            fsi::matrix smatEmbedded;
            fsi::matrix qmatEmbedded;
            fsi::vector dsdc;

            // Store function in memory in case the source term is requested
            // by the solver
            bool corrector;
            int stageIndex;
            fsi::matrix Sj;
            bool convergence;
            int timeIndex;
            const int minSweeps;
            const int maxSweeps;

            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
            std::shared_ptr<sdc::DataStorage> data;
    };
}

#endif
