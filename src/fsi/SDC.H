
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

namespace sdc
{
    class SDC
    {
public:

        SDC(
            std::shared_ptr<SDCSolver> solver,
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
            std::string rule,
            int nbNodes,
            double tol
            );

        SDC(
            std::string rule,
            int nbNodes,
            double tol
            );

        ~SDC();

        void run();

        void solveTimeStep( const double t0 );

        void getSourceTerm(
            const bool corrector,
            const int k,
            const double deltaT,
            Eigen::VectorXd & rhs,
            Eigen::VectorXd & qold
            );

        void setFunction(
            const int k,
            const Eigen::VectorXd & f,
            const Eigen::VectorXd & result
            );

        void setOldSolution( const Eigen::VectorXd & result );

        void outputResidual( std::string name );

        bool isConverged();

        std::shared_ptr<SDCSolver> solver;
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper;

        int nbNodes;
        int N;
        int k;

        double dt;
        double tol;

        Eigen::VectorXd nodes;
        Eigen::MatrixXd smat;
        Eigen::MatrixXd qmat;
        Eigen::VectorXd nodesEmbedded;
        Eigen::MatrixXd smatEmbedded;
        Eigen::MatrixXd qmatEmbedded;
        Eigen::VectorXd dsdc;

private:

        void computeResidual(
            const Eigen::MatrixXd & qmat,
            const Eigen::MatrixXd & F,
            const double dt,
            Eigen::MatrixXd & qj
            );

        // Store function in memory in case the source term is requested
        // by the solver
        bool corrector;
        int stageIndex;
        double stageDeltaT;
        Eigen::MatrixXd F;
        Eigen::MatrixXd Fold;
        Eigen::MatrixXd Sj;
        Eigen::MatrixXd solStages;
        bool convergence;
    };
}

#endif
