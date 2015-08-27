
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Oscillator_H
#define Oscillator_H

#include "ESDIRK.H"
#include "SDC.H"
#include <cmath>

using namespace sdc;

class Oscillator : public SDCSolver
{
public:

    Oscillator(
        int nbTimeSteps,
        scalar dt,
        fsi::vector q0,
        scalar amplitude,
        scalar frequency,
        scalar m,
        scalar k
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

    fsi::vector sol;

private:

    int nbTimeSteps;
    scalar dt;
    fsi::vector q0;
    scalar amplitude;
    scalar frequency;
    scalar m;
    scalar k;
    scalar t;
};

#endif
