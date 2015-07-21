
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Cos_H
#define Cos_H

#include "ESDIRK.H"
#include "SDC.H"

using namespace sdc;

class Cos : public SDCSolver
{
public:

    Cos(
        int nbTimeSteps,
        double dt,
        double endTime,
        double amplitude,
        double frequency
        );

    virtual void evaluateFunction(
        const int k,
        const Eigen::VectorXd & q,
        const double t,
        Eigen::VectorXd & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual double getScalingFactor();

    virtual void getSolution( Eigen::VectorXd & solution );

    virtual void setSolution(
        const Eigen::VectorXd & solution,
        const Eigen::VectorXd & f
        );

    virtual double getEndTime();

    virtual double getTimeStep();

    virtual void nextTimeStep();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void setNumberOfImplicitStages( int k );

    virtual void implicitSolve(
        bool corrector,
        const int k,
        const int kold,
        const double t,
        const double dt,
        const Eigen::VectorXd & qold,
        const Eigen::VectorXd & rhs,
        Eigen::VectorXd & f,
        Eigen::VectorXd & result
        );

    double f;
    double sol;

private:

    int nbTimeSteps;
    double dt;
    double endTime;
    double amplitude;
    double frequency;
    double t;
};

#endif
