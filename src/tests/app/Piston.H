
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Piston_H
#define Piston_H

#include "ESDIRK.H"

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
        double dt,
        double q0,
        double qdot0,
        double As,
        double Ac,
        double omega
        );

    ~Piston();

    double referenceSolution( double t );

    void run();

    virtual void finalizeTimeStep(){}

    virtual void getSolution( Eigen::VectorXd & solution );

    virtual void initTimeStep(){}

    virtual bool isRunning();

    virtual void implicitSolve(
        bool corrector,
        const int k,
        const double t,
        const double dt,
        const Eigen::VectorXd & qold,
        const Eigen::VectorXd & rhs,
        Eigen::VectorXd & f,
        Eigen::VectorXd & result
        );

    virtual int getDOF();

    virtual double getTimeStep();

    virtual void evaluateFunction(
        const int k,
        const Eigen::VectorXd & q,
        const double t,
        Eigen::VectorXd & f
        );

    virtual void setDeltaT( double dt ){}

    virtual void setNumberOfStages( int k ){}

    virtual void nextTimeStep();

    virtual int getNbCells()
    {
        return 2;
    }

    virtual double getScalingFactor()
    {
        return 1;
    }

    int nbTimeSteps;
    double dt;
    double q0;
    double c1;
    double c2;
    double As;
    double Ac;
    double omega;
    int N;
    double q;
    double qdot;
    double t;
    double timeIndex;
};

#endif
