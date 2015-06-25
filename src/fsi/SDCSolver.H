
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SDCSolver_H
#define SDCSolver_H

#include <Eigen/Dense>

namespace sdc
{
    class SDCSolver
    {
public:

        virtual ~SDCSolver(){}

        virtual void evaluateFunction(
            const int k,
            const Eigen::VectorXd & q,
            const double t,
            Eigen::VectorXd & f
            ) = 0;

        virtual void finalizeTimeStep() = 0;

        virtual int getDOF() = 0;

        virtual double getScalingFactor() = 0;

        virtual void getSolution( Eigen::VectorXd & solution ) = 0;

        virtual void setSolution(
            const Eigen::VectorXd & solution,
            const Eigen::VectorXd & f
            ) = 0;

        virtual double getEndTime() = 0;

        virtual double getTimeStep() = 0;

        virtual void nextTimeStep() = 0;

        virtual void initTimeStep() = 0;

        virtual bool isRunning() = 0;

        virtual void setDeltaT( double dt ) = 0;

        virtual void setNumberOfStages( int k ) = 0;

        virtual void implicitSolve(
            bool corrector,
            const int k,
            const double t,
            const double dt,
            const Eigen::VectorXd & qold,
            const Eigen::VectorXd & rhs,
            Eigen::VectorXd & f,
            Eigen::VectorXd & result
            ) = 0;
    };
}

#endif
