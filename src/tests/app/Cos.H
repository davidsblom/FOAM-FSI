
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
        scalar dt,
        scalar endTime,
        scalar amplitude,
        scalar frequency
        );

    virtual void evaluateFunction(
        const int k,
        const fsi::vector & q,
        const scalar t,
        fsi::vector & f
        );

    virtual void finalizeTimeStep();

    virtual int getDOF();

    virtual scalar getScalingFactor();

    virtual void getSolution( fsi::vector & solution );

    virtual void setSolution(
        const fsi::vector & solution,
        const fsi::vector & f
        );

    virtual scalar getEndTime();

    virtual scalar getTimeStep();

    virtual void nextTimeStep();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void setNumberOfImplicitStages( int k );

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
        );

    scalar f;
    scalar sol;

private:

    int nbTimeSteps;
    scalar dt;
    scalar endTime;
    scalar amplitude;
    scalar frequency;
    scalar t;
};

#endif
