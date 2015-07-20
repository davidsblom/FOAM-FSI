
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

namespace sdc
{
    class ESDIRK
    {
public:

        ESDIRK(
            std::shared_ptr<SDCSolver> solver,
            std::string method,
            std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper
            );

        ESDIRK( std::string method );

        ~ESDIRK();

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

        int getNbStages();

        bool isStageImplicit( double Akk );

        void setOldSolution( const Eigen::VectorXd & result );

        std::shared_ptr<SDCSolver> solver;
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper;
        double dt;
        int nbStages;
        Eigen::MatrixXd A;
        Eigen::VectorXd B;
        Eigen::VectorXd C;
        Eigen::VectorXd Bhat;
        int N;

private:

        void initializeButcherTableau( std::string method );

        // Store function in memory in case the source term is requested
        // by the solver. This functionality is used by the moving wall
        // velocity boundary condition.
        int stageIndex;
        Eigen::MatrixXd F;
        Eigen::MatrixXd solStages;
        Eigen::VectorXd qold;
    };
}

#endif
