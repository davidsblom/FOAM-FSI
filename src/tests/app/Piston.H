
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Piston_H
#define Piston_H

#include "ESDIRK.H"
#include "SDC.H"

using namespace sdc;

/*
 * Implements the piston problem for which a exact solution exists.
 * With the exact solution, the convergence orders can be confirmed
 * exactly for this test case.
 *
 * Reference: A. H. van Zuijlen, Efficient Higher Order Time
 * Integration of Partitioned Systems, PhD thesis, page 22.
 *
 * Notes Gaussian quadrature:
 * http://undergraduate.csse.uwa.edu.au/units/GENG2140/Mech2012/Gaussian_Quadrature.pdf
 * http://en.wikipedia.org/wiki/Gaussian_quadrature
 * http://en.wikipedia.org/wiki/Legendre_polynomials
 */

class Piston : public SDCSolver
{
public:

    Piston(
        int nbTimeSteps,
        scalar dt,
        scalar q0,
        scalar qdot0,
        scalar As,
        scalar Ac,
        scalar omega
        );

    Piston(
        int nbTimeSteps,
        scalar dt,
        scalar q0,
        scalar qdot0,
        scalar As,
        scalar Ac,
        scalar omega,
        std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme,
        int k
        );

    ~Piston();

    scalar referenceSolution( scalar t );

    void run();

    virtual void finalizeTimeStep(){}

    virtual void getSolution(
        fsi::vector & solution,
        fsi::vector & f
        );

    virtual void setSolution(
        const fsi::vector & solution,
        const fsi::vector & f
        );

    virtual scalar getEndTime();

    virtual void initTimeStep(){}

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

    virtual int getDOF();

    virtual scalar getTimeStep();

    virtual void evaluateFunction(
        const int k,
        const fsi::vector & q,
        const scalar t,
        fsi::vector & f
        );

    virtual void setNumberOfImplicitStages( int k );

    virtual void nextTimeStep();

    virtual scalar getScalingFactor()
    {
        return 1;
    }

    virtual void getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        );

    int nbTimeSteps;
    scalar dt;
    scalar q0;
    scalar c1;
    scalar c2;
    scalar As;
    scalar Ac;
    scalar omega;
    int N;
    scalar q;
    scalar qdot;
    scalar t;
    scalar timeIndex;
    scalar endTime;
    int k;

    std::deque<fsi::vector> solStages;
    std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme;
};

#endif
